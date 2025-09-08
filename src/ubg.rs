use crate::rearrangement::*;
use crate::gfa::{LEN_PREFIX,get_or_set_node_id,parse_marker};
use crate::util::{reverse_map, SAFE_GFA_OVERLAP};
use std::collections::{HashMap, HashSet};
use csv::{ReaderBuilder};
use std::io::{self, BufRead, BufReader};
use std::fs::File;


#[derive(Debug,Clone)]
pub struct UBG {
    pub node_sizes : HashMap<Marker,usize>,
    pub adjacencies : HashMap<Extremity,HashSet<Extremity>>,
    pub node_ids : HashMap<String,Marker>
}



fn remove_node_from_adj(adjacencies : &mut HashMap<Extremity,HashSet<Extremity>>, xtr : Extremity) {
    let neighbors=adjacencies.get(&xtr).cloned().expect("aa");
    for x in neighbors.iter() {
        let v: &mut HashSet<Extremity> =  adjacencies.get_mut(&x).expect("A");
        v.remove(&xtr);
    }
    adjacencies.remove(&xtr);
}

impl RearrangementGraph for UBG  {
    fn degree(&self,n:Extremity) -> Option<usize> {
        let neighb = self.adjacencies.get(&n)?;
        if ! neighb.contains(&n) {
            return Some(neighb.len())
        }
        Some(neighb.len()+1)
    }

    fn adj_neighbors(&self,n:Extremity) -> Option<impl Iterator<Item=Extremity>> {
        Some(self.adjacencies.get(&n)?.iter().copied())
    }
    
    fn node_size(&self,n:Marker) -> Option<usize> {
        Some(*self.node_sizes.get(&n)?)
    }
    fn markers(&self) -> impl Iterator<Item=Marker> {
        self.node_ids.values().into_iter().copied()
    }

    fn trim_singlethread(&mut self, threshold : usize) {
        let mut to_remove = HashSet::new();
        for (node,sz) in self.node_sizes.iter() {
            if *sz < threshold {
                //println!("Removing {} , i.e. {} (hd) {} (tl)",node,head(*node),tail(*node));
                let hd = head(*node);
                let tl = tail(*node);
                let nh = self.adjacencies.get(&hd).expect("Assertion violated: Marker extremity not in adjacencies.").clone();
                let nt = self.adjacencies.get(&tl).expect("Assertion violated: Marker extremity not in adjacencies.").clone();
                //add adjacencies between the neighboring markers
                for x in nh.iter() {
                    for y in nt.iter() {
                        self.adjacencies.get_mut(&x).expect("!").insert(*y);
                        self.adjacencies.get_mut(&y).expect("!").insert(*x);
                    }
                }
                if nh.contains(&hd) {
                    for x in nt.iter() {
                        for y in nt.iter() {
                            self.adjacencies.get_mut(&x).expect("!").insert(*y);
                            self.adjacencies.get_mut(&y).expect("!").insert(*x);
                        }
                    }
                }
                if nt.contains(&tl) {
                    for x in nh.iter() {
                        for y in nh.iter() {
                            self.adjacencies.get_mut(&x).expect("!").insert(*y);
                            self.adjacencies.get_mut(&y).expect("!").insert(*x);
                        }
                    }
                }
                remove_node_from_adj(&mut self.adjacencies, hd);
                remove_node_from_adj(&mut self.adjacencies, tl);
                to_remove.insert(*node);   
            }
        }

        for node in &to_remove {
            self.node_sizes.remove(node);
        }
        let mut xd = Vec::new();
        for (x,y) in &self.node_ids{
            if to_remove.contains(y) {
                xd.push(x.clone());
            }
        }

        for x in xd {
            self.node_ids.remove(&x);
        }
        //remove any accidentally created self loops of the telomere
        if let Some(tladj) = self.adjacencies.get_mut(&TELOMERE) {
            tladj.remove(&TELOMERE);
        }

        //remove the telomere if it's empty
        if self.degree(TELOMERE).unwrap_or(0) == 0 {
            self.adjacencies.remove(&TELOMERE);
        }
    }
    
    fn trim_multithread(&mut self, min_size : usize,n_threads : usize) {
        if n_threads <= 1 {
            self.trim_singlethread(min_size);
        } else {
             panic!("Not implemented!");
        }
       
    }


    fn iter_adjacencies(&self) -> impl Iterator<Item=Adjacency> {
        self.adjacencies.iter().flat_map(|(x,neighbors)| {
            neighbors.iter().filter_map(|y| {
                if *y >= *x {
                    Some((*x,*y))
                } else {
                    None
                }
            })
        })
    }

    fn extremities(&self) -> impl Iterator<Item=Extremity> {
        self.adjacencies.keys().into_iter().copied()
    }

    fn num_extremities(&self) -> usize {
        self.adjacencies.len()
    }

    fn num_markers(&self) -> usize {
        self.node_ids.len()
    }

    fn from_hash_maps(sizes : HashMap<Marker,usize>, adj :HashMap<Extremity,HashSet<Extremity>> ,  nids : HashMap<String,Marker>) -> Self {
        return UBG { node_sizes:sizes, adjacencies: adj, node_ids: nids }
    }

    fn from_gfa(path: &str,ignore_overlap : bool) -> io::Result<UBG>{
    if !ignore_overlap {
        panic!("Not implemented.");
    }
    eprintln!("Read gfa.");
    let mut node_sizes = HashMap::new();
    let mut adjacencies = HashMap::new(); 
    let mut node_ids: HashMap<String, Marker>   = HashMap::new();
    let mut rdr = ReaderBuilder::new().has_headers(false).delimiter(b'\t').flexible(true).from_path(path)?;
    let mut curr_id = 1;
    let mut i :u32 = 0;
    let mut n_edges :usize = 0;
    for res in rdr.records() {
        let x = res?;
        i+=1;
        if i%1000000==0 {
            eprintln!("Read {} lines.",i);
            eprintln!("{} nodes and {} edges in graph.",node_sizes.len(),n_edges);
        }
  
        let entrytype= match x.get(0) {
            None => continue,
            Some(y) => y
        };
        
        if entrytype=="S" {
            let seg_name = x.get(1);
            let seg_str = x.get(2);
            let seg_name= match seg_name {
                None => return Err(io::Error::new(io::ErrorKind::Other,"Empty segment label")),
                Some(y) => y.to_string()
            };
            let mut seg_len = match seg_str {
                None => 0,
                Some(y) => y.len()
            };
            for entry in x.iter().skip(3) {
                if entry.starts_with(&LEN_PREFIX) {
                    let lenstr = &entry[LEN_PREFIX.len()..];
                    seg_len = lenstr.parse().expect(&format!("Invalid length: {}",lenstr));
                }
            }
            let n_id;
            (curr_id,n_id) = get_or_set_node_id(&mut node_ids, curr_id, seg_name);
            if !node_sizes.contains_key(&n_id) {
                node_sizes.insert(n_id, seg_len);
            }
            
            
        } else if entrytype == "L" {
            let (sega,segb) = match  (x.get(1),x.get(3)) {
                (Some(a),Some(b)) => (a,b),
                (_,_) => return Err(io::Error::new(io::ErrorKind::Other,format!("Malformed link: {:?} in line {}",x,i)))
            };
            let aid;
            let bid;
            (curr_id,aid) = get_or_set_node_id(&mut node_ids, curr_id, sega.to_string());
            (curr_id,bid) = get_or_set_node_id(&mut node_ids, curr_id, segb.to_string());
            let (axtr,bxtr)= match (x.get(2),x.get(4)) {
                (Some("+"),Some("+")) => (head(aid),tail(bid)),
                (Some("+"),Some("-")) => (head(aid),head(bid)),
                (Some("-"),Some("+")) => (tail(aid),tail(bid)),
                (Some("-"),Some("-")) => (tail(aid),head(bid)),
                (_,_) => return Err(io::Error::new(io::ErrorKind::Other,format!("Malformed link: {:?} in line {}",x,i)))
            };
            if !adjacencies.contains_key(&axtr) {
                adjacencies.insert(axtr,HashSet::new());
            }
            if !adjacencies.contains_key(&bxtr) {
                adjacencies.insert(bxtr,HashSet::new());
            }
            adjacencies.get_mut(&axtr).expect("Horror").insert(bxtr);
            adjacencies.get_mut(&bxtr).expect("Horror").insert(axtr);
            n_edges+=1;
        }

    }
    println!("GFA file parsing done.");
    println!("Read {} lines.",i);
    println!("{} nodes and {} edges in graph.",node_sizes.len(),n_edges);
    
    Ok(UBG {
        node_sizes,
        adjacencies,
        node_ids
    })
}

fn fill_telomeres (&mut self) {
    for (node,_) in self.node_sizes.clone().iter() {
        check_add_tel(self, head(*node));
        check_add_tel(self, tail(*node));
    }
}

fn from_unimog(path : &str) -> io::Result<UBG> {
    let node_sizes = HashMap::new();
    let mut adjacencies :HashMap<Extremity, HashSet<Extremity>> = HashMap::new();
    let mut node_ids: HashMap<String, Marker>   = HashMap::new();
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut curr_id = 1;
    for line in reader.lines() {
        let line  = &line?;
        if line.starts_with(&">") {
            continue;
        }
        let mut line = &line[..];
        let mut first = None;
        let mut last = None;
        while let Some(j) =  line.find(" ") {
            let m = &line[0..j];
            
            eprintln!("{}",m);
            let is_forward;
            let marker;
            (curr_id,is_forward,marker) = parse_marker(&mut node_ids, m, curr_id);
            //eprintln!("is_forward {}",is_forward);
            if first.is_none() {
                first=Some((is_forward,marker));
            }
            if let Some(last) = last {
                let (xta,xtb) = to_adjacency(last,(is_forward,marker));
                adjacencies.entry(xta).or_insert(HashSet::new()).insert(xtb);
                adjacencies.entry(xtb).or_insert(HashSet::new()).insert(xta);
            }
            line = &line[j+1..];
            last = Some((is_forward,marker));
            if line.starts_with(")") || line.starts_with("|") {
                break;
            }
        }
        if let Some((ifa,ma)) = last {
            let (ifb,mb) = first.expect("If a last marker exists, there must be a first marker");
            let (xta,xtb) = to_adjacency((ifa,ma),(ifb,mb));
            if line.starts_with(")") {
                adjacencies.entry(xta).or_insert(HashSet::new()).insert(xtb);
                adjacencies.entry(xtb).or_insert(HashSet::new()).insert(xta);
            } else if line.starts_with("|") {
                adjacencies.entry(TELOMERE).or_insert(HashSet::new()).insert(xtb);
                adjacencies.entry(TELOMERE).or_insert(HashSet::new()).insert(xta);
                adjacencies.entry(xta).or_insert(HashSet::new()).insert(TELOMERE);
                adjacencies.entry(xtb).or_insert(HashSet::new()).insert(TELOMERE);
            } else {
                return Err(io::Error::new(io::ErrorKind::Other,"Invalid chromosome end."));
            }
        }

    }
    Ok(UBG {
        node_sizes,
        adjacencies,
        node_ids
    })

   
}

fn name_to_marker(&self,name : &str) -> Option<Marker> {
    self.node_ids.get(name).copied()
}


fn marker_names(&self) -> HashMap<Marker,String> {
    reverse_map(&self.node_ids)
}

fn overlap(&self,_:Extremity,_:Extremity) -> usize {
    return 0;
}
}


fn check_add_tel(ubg : &mut UBG, n : Extremity) {
    match ubg.adjacencies.get(&n) {
        None => ubg.adjacencies.insert(n, HashSet::new()),
        Some(_) => None
    };
    let neighbors= ubg.adjacencies.get_mut(&n).expect("Just inserted value disappeared");
    if neighbors.is_empty() {
        neighbors.insert(TELOMERE);
        match ubg.adjacencies.get(&TELOMERE) {
            None => ubg.adjacencies.insert(TELOMERE, HashSet::new()),
            Some(_) => None
        };
        ubg.adjacencies.get_mut(&TELOMERE).expect("test").insert(n);
    }
    
}