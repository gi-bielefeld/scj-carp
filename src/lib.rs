use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::hash::Hash;
use std::io::{self, BufRead, BufReader};
use csv::{ReaderBuilder,Reader};
use std::cmp::Ordering;
use std::thread;
use flate2::read::MultiGzDecoder;
use itertools::Itertools;


const ONE_MILLION :usize = 1000000;
const LEN_PREFIX  : &str = "LN:i:";
const TELOMERE : Extremity = 0;

pub const CARP_LOGO : &str = 
"   oo.              You're running SCJ CARP
 o00O                         Version 0.0.1
  O00O                                     
  O00.                                     
  O0O          Feel free to report bugs at:
  o0o      github.com/gi-bielefeld/scj-carp
   0o                                      
   O0O                                     
   O00O.                                   
   o0000Oo                                 
   .0000000Oo.                             
    o000000000Oo.                          
     O00000000000Oo.                       
      O0000000000000Oo.                    
       o0000000000000000Oo.                
        o000000000000000000Oo.             
         .O0000000000000000000o.           
           o00000000000000000000Oo.        
            .O000000000000Oo....o00o.      
           .OO000000000000o......O000.     
          .0000000000000000OoooO000000.    
          .O000000OO0000000000000000000.   
            ..ooo.   .O00000000000000000.  
                  .oO00OOO000000000000O... 
               .oOOo..    ..ooooOOO0o.     
              ..                   .       
";


pub const CARP_VERSION : &str = "0.0.1";

#[derive(Debug,Clone)]
pub struct UBG {
    pub node_sizes : HashMap<Marker,usize>,
    pub adjacencies : HashMap<Extremity,HashSet<Extremity>>,
    pub node_ids : HashMap<String,Marker>
}


pub trait RearrangementGraph
where Self : Sized,
Self: Send,
Self: Sync
{
    fn degree(&self,n:Extremity) -> Option<usize>;
    fn adj_neighbors(&self,n:Extremity) -> Option<impl Iterator<Item=Extremity>>;
    fn node_size(&self,n:Marker) -> Option<usize>;
    fn trim_singlethread(&mut self, min_size : usize);
    fn trim_multithread(&mut self, min_size : usize,n_threads : usize);
    fn markers(&self) -> impl Iterator<Item=Marker>;
    fn iter_adjacencies(&self) -> impl Iterator<Item=Adjacency>;
    fn extremities(&self) -> impl Iterator<Item=Extremity>;
    fn num_markers(&self) -> usize;
    fn num_extremities(&self) -> usize;
    fn from_hash_maps(sizes : HashMap<Marker,usize>, adj :HashMap<Extremity,HashSet<Extremity>> ,  nids : HashMap<String,Marker>) -> Self;
    fn from_gfa(path: &str) -> io::Result<Self>;
    fn fill_telomeres(&mut self);
    fn from_unimog(path : &str) -> io::Result<Self>;
    fn name_to_marker(&self,name : &str) -> Option<Marker>;
    fn marker_names(&self) -> HashMap<Marker,String>;
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

    fn from_gfa(path: &str) -> io::Result<UBG>{
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
}


#[derive(Debug,Clone)]
pub struct MBG {
    node_sizes : Vec<usize>,
    adjacencies : Vec<Vec<Extremity>>,
    node_ids : HashMap<String,Marker>,
    masked_markers : HashSet<Marker>,
}


impl MBG {
    #[inline(always)]
    fn remove_marker(&mut self, m : Marker) {
        if self.masked_markers.contains(&m) {
            return
        }
        let tailnb : HashSet<Extremity> = self.adjacencies.get(tail(m)).unwrap().iter().copied().collect();
        let headnb : HashSet<Extremity>  = self.adjacencies.get(head(m)).unwrap().iter().copied().collect();
        self.adjacencies[head(m)]= Vec::new();
        self.adjacencies[tail(m)] = Vec::new();
        //It's a bit annoying going with the HashMap intermediary
        for x in &tailnb {
            if *x==head(m) || *x==tail(m) {
                continue;
            }
            let mut xneighbors : HashSet<Extremity> = self.adjacencies.get(*x).unwrap().iter().copied().collect();
            xneighbors.remove(&tail(m));
            for y in &headnb {
                if *y != head(m) && *y != tail(m) {
                    xneighbors.insert(*y);
                }
            }
            if headnb.contains(&head(m)) {
                for y in &tailnb {
                    if *y != head(m) && *y != tail(m) {
                        xneighbors.insert(*y);
                    }
                }
            }
            //prevent 0,0 adjacency
            if *x==TELOMERE && xneighbors.contains(&TELOMERE) {
                xneighbors.remove(&TELOMERE);
            }
            self.adjacencies[*x] = xneighbors.iter().copied().collect();
            //for degree purposes add x to itself twice
            if xneighbors.contains(x) {
                self.adjacencies[*x].push(*x);
            }
        }

        for x in &headnb {
            if *x==head(m) || *x==tail(m) {
                continue;
            }
            let mut xneighbors : HashSet<Extremity> = self.adjacencies.get(*x).unwrap().iter().copied().collect();
            xneighbors.remove(&head(m));
            for y in &tailnb {
                if *y != head(m) && *y != tail(m) {
                    xneighbors.insert(*y);
                }
            }
            if tailnb.contains(&tail(m)) {
                for y in &headnb {
                    if *y != head(m) && *y != tail(m) {
                        xneighbors.insert(*y);
                    }
                } 
            }
             //prevent 0,0 adjacency
            if *x==TELOMERE && xneighbors.contains(&TELOMERE) {
                xneighbors.remove(&TELOMERE);
            }
            self.adjacencies[*x]=xneighbors.iter().copied().collect();
            //for degree purposes add x to itself twice
            if xneighbors.contains(x) {
                self.adjacencies[*x].push(*x);
            }
        }

        


        self.masked_markers.insert(m); 
    }

    fn fill_up_vec<T>(adj : &mut Vec<T>, to : usize)
    where T: std::default::Default
     {
        for _ in adj.len()..(to+1) {
            adj.push(T::default());
        }
    }
    
    
    fn find_solid_neighbors(&self, start : Extremity, size_threshold : usize) -> HashSet<Extremity> {
        let mut stack = Vec::new();
        let mut visited = HashSet::new();
        let mut solid_neighbors = HashSet::new();
        stack.push(start);
        while stack.len() > 0 {
            let x = stack.pop().unwrap();
            visited.insert(x);
            for y in self.adj_neighbors(x).unwrap() {
                let z = other(y);

                    if self.node_size(marker(y)).unwrap_or(0) >= size_threshold || z==TELOMERE {
                        solid_neighbors.insert(y);
                    } else if !visited.contains(&z) {
                        //eprintln!("Skipping over marker {}, size {:?}",marker(y),self.node_size(marker(y)));
                        stack.push(z);
                    }
                
            }
        }
        solid_neighbors
    }


    fn gfa_from_any_reader<R>(rdr : &mut Reader<R>) -> Result<MBG, io::Error>
    where 
        R : std::io::Read
    {
        eprintln!("Read gfa.");
        let mut node_sizes = Vec::new();
        let mut adjacencies : Vec<Vec<Extremity>> = Vec::new(); 
        let mut node_ids: HashMap<String, Marker>   = HashMap::new();
        let mut curr_id = 1;
        let mut i :usize = 0;
        let mut n_edges :usize = 0;
        let mut telomeres : HashSet<(String,bool)> = HashSet::new();
        let mut seen_edges = HashSet::new();
        for res in rdr.records() {
            let x = res?;
            i+=1;
            if i%ONE_MILLION==0 {
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
                if n_id >= node_sizes.len() {
                    Self::fill_up_vec(&mut node_sizes, n_id);
                }
                node_sizes[n_id] = seg_len;
                
                
                
            } else if entrytype == "L" {
                let (sega,segb) = match  (x.get(1),x.get(3)) {
                    (Some(a),Some(b)) => (a,b),
                    (_,_) => return Err(io::Error::new(io::ErrorKind::Other,"Malformed link."))
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
                    (_,_) => return Err(io::Error::new(io::ErrorKind::Other,"Malformed link."))
                };
                if adjacencies.len() <= axtr {
                    Self::fill_up_vec(&mut adjacencies, axtr);
                }
                if adjacencies.len() <= bxtr {
                    Self::fill_up_vec(&mut adjacencies, bxtr);
                }
                //adjacencies.get_mut(&axtr).expect("Horror").insert(bxtr);
                //adjacencies.get_mut(&bxtr).expect("Horror").insert(axtr);
                let cane = canonicize((axtr,bxtr));
                if !seen_edges.contains(&cane) {
                    adjacencies.get_mut(axtr).unwrap().push(bxtr);
                    adjacencies.get_mut(bxtr).unwrap().push(axtr);
                    n_edges+=1;
                    seen_edges.insert(cane);
                }
                
            } else if entrytype == "P" {
                let pname = x.get(1).expect("Path does not have a name identifier.");
                let mut path = x.get(2).expect(&format!("Path '{pname}' missing mandatory gfa field 3.")).split(|x : char| {x==',' || x==';'});
                let fst = path.nth(0);
                //let lst = path.last();
                let lst = match path.last() {
                    Some(a) => Some(a),
                    None => fst
                };
                let parse_pend = |x : &str,is_path_end : bool| {
                    let x = x.strip_suffix("\n").unwrap_or(x);
                    assert!(x.ends_with("+") || x.ends_with("-"));
                    let mut xp = x.to_owned();
                    xp.pop();
                    let xtr_is_tail = is_path_end == x.ends_with("-");
                    (xp,xtr_is_tail)};
                match (fst,lst) {
                    (Some(f),Some(l)) => {
                    telomeres.insert(parse_pend(f,false));
                    telomeres.insert(parse_pend(l,true));
                    },
                    (_,_) => continue
                }
            } else if entrytype == "W" {
                let wlk = x.get(6).expect("Walk line without walk");
                let pat = |x : char| {x=='>' || x=='<'};
                let start = wlk.find(pat);
                let end = wlk.rfind(pat);
                let mut wlki = wlk[1..].split(pat);
                let fst = wlki.nth(0);
                let lst = match wlki.last() {
                    Some(a) => Some(a),
                    None => fst
                };
                match (fst,lst) {
                    (Some(f),Some(l)) => {
                    let e_is_tail = wlk.as_bytes()[end.unwrap()] as char == '<';
                    let s_is_tail = wlk.as_bytes()[start.unwrap()] as char == '>';
                    telomeres.insert((f.to_owned(),s_is_tail));
                    telomeres.insert((l.to_owned(),e_is_tail));
                    },
                    (_,_) => continue
                }
            }

        }
        eprintln!("Filling in {} telomeres observed in paths",telomeres.len());
        for (mrk,xtr_is_tail) in telomeres {
            let m = *node_ids.get(&mrk).expect(&format!("Segment {mrk} occurs in a path, but not as a segment entry."));
            let xtr = if xtr_is_tail {
                tail(m)
            } else  {
                head(m)
            };
            if adjacencies.len() <= xtr {
                    Self::fill_up_vec(&mut adjacencies, xtr);
            }
            adjacencies.get_mut(xtr).unwrap().push(TELOMERE);
            adjacencies.get_mut(TELOMERE).unwrap().push(xtr);
        }
        //eprintln!("{node_sizes:?}");
        
        Ok(MBG { node_sizes: node_sizes, adjacencies: adjacencies, node_ids: node_ids, masked_markers: HashSet::from([TELOMERE]) })
    }
    
    fn identify_removal_nodes_in_range(&self, min_size: usize,from : Marker, to:Marker) -> Vec<Marker> {
        let mut to_remove = Vec::new();
        if from >=to {
            return  to_remove;
        }
        for (offset,size) in self.node_sizes[from..to].iter().enumerate() {
            let m = offset+from;
            if *size < min_size {
                to_remove.push(m);
            }
        }
        to_remove
    }

    fn identify_removal_nodes(&self, min_size: usize) -> Vec<Marker> {
        self.identify_removal_nodes_in_range(min_size, 1, self.num_markers()+1)
    }



    fn identify_removal_nodes_mthread(&self,min_size: usize,n_threads : usize) -> Vec<Marker> {
        let mut to_remove = Vec::new();
        thread::scope(|scope| {
            let nmarkers =  self.num_markers();
            let mut handles = Vec::new();
            let slice_size =nmarkers/n_threads +1;
            for i in 0..n_threads {
                let g = &self;
                let lb = (slice_size*i).min(nmarkers);
                let rb = (slice_size*(i+1)).min(nmarkers);
                eprintln!("Spawning thread {i} processing markers with index {lb} to {rb} (total {nmarkers})");
                let x =  scope.spawn(move || g.identify_removal_nodes_in_range(min_size, lb, rb));
                handles.push(x);
            }
            eprintln!("Joining results.");
            for x in handles {
                to_remove.extend(x.join().unwrap());
            }
        
    });
        to_remove
    }
    
    fn remove_all(&mut self, to_remove: &Vec<Marker>) {
        let tsize = to_remove.len();
        let about_one_percent = (tsize/100).max(1).min(ONE_MILLION/10);
        for (i,m) in to_remove.iter().copied().enumerate() {
            if (i+1)%about_one_percent==0{
                eprintln!("Removal processed {i}/{tsize} markers.");
            }
            self.remove_marker(m);

        }
    }

    pub fn trim_any(&mut self,min_size: usize, n_threads : usize) {
        if n_threads == 1 {
            eprintln!("Trimming singlethreaded.");
            self.trim_singlethread(min_size);
            return;
        }
        let rm = self.identify_removal_nodes_mthread(min_size, n_threads);
        if rm.len() >= self.num_markers()/10 {
            eprintln!("More then 10% of nodes are scheduled for removal. Trimming singlethreaded.");
            self.remove_all(&rm);
            return;
        }
        self.trim_multithread(min_size, n_threads);
    }

}





fn has_deleted_neighbor(graph : &MBG, extremity : Extremity, size_threshold : usize) -> bool {
    for x in graph.adj_neighbors(extremity).unwrap() {
        if graph.node_size(marker(x)).unwrap_or(0) < size_threshold && x != TELOMERE {
            return true;
        }
    }
    false
}

fn __trim_vertices(graph : &MBG, from : Extremity, to : Extremity, size_threshold : usize) -> (HashSet<Marker>,Vec<(Extremity,Vec<Extremity>)>) {
        let mut masked = HashSet::new();
        let mut adjacencies = Vec::new();
        let nnodes = to -from;
        //eprintln!("Size thresh {size_threshold}");
        for i in from..to {
            let m = marker(i);
            if i == 1 {
                continue;
            }
            if graph.masked_markers.contains(&m) && i != TELOMERE {
                continue;
            }
            if graph.node_size(m).unwrap_or(0) < size_threshold && m != TELOMERE {
                masked.insert(m);
                //eprintln!("Deleting marker {m}...");
            } else if has_deleted_neighbor(graph, i, size_threshold) {
                let mut sn = graph.find_solid_neighbors(i, size_threshold);
                //prevent adjacency 0,0 from appearing
                if i==TELOMERE && sn.contains(&TELOMERE) {
                    sn.remove(&TELOMERE);
                }
                let mut locadj : Vec<Extremity> = sn.iter().copied().collect();
                if sn.contains(&i) {
                    //degree purposes
                    locadj.push(i);
                }
                adjacencies.push((i,locadj));
            }
            if (i - from)%ONE_MILLION == 0 && i > from {
                let proc_nodes = i-from;
                let percentage = proc_nodes as f64 / nnodes as f64 * 100.0;
                eprintln!("Processed {proc_nodes} extremities in range [{from},{to}[ ({percentage:.02}%)");
            }
        }
        //eprintln!("{from} {to}: {adjacencies:?}");
        (masked,adjacencies)
    }

impl RearrangementGraph for MBG {
    fn degree(&self,n:Extremity) -> Option<usize> {
        if self.masked_markers.contains(&marker(n)) && n!=TELOMERE {
            return None;
        }
        self.adjacencies.get(n).map(|x| x.len())
    }

    fn adj_neighbors(&self,n:Extremity) -> Option<impl Iterator<Item=Extremity>> {
        if self.masked_markers.contains(&marker(n)) && n!=TELOMERE {
            return None;
        }
        self.adjacencies.get(n).map(|x| x.iter().copied())
    }

    fn markers(&self) -> impl Iterator<Item=Marker> {
        self.node_ids.values().filter_map(|x| {if self.masked_markers.contains(x) {
            None
        } else {
            Some(*x)
        }
        })
    }

    fn extremities(&self) -> impl Iterator<Item=Extremity> {
        self.markers().flat_map(|m| [tail(m),head(m)]).filter(|x| !self.masked_markers.contains(&marker(*x))).chain([TELOMERE])
    }

    fn iter_adjacencies(&self) -> impl Iterator<Item=Adjacency> {
        self.adjacencies.iter().enumerate().flat_map(|(xtr,neighb)|{
         neighb.iter().filter_map(move |ytr| {
        if xtr <= *ytr {
            Some((xtr,*ytr))
        } else {
            None
            }})    
        }).unique()
    }

    fn node_size(&self,n:Marker) -> Option<usize> {
        if n == TELOMERE || self.masked_markers.contains(&n) {
            return None;
        }
        self.node_sizes.get(n).copied()
    }

    fn num_markers(&self) -> usize {
        //+1 because of telomere node, which does not have an id, but occurs in masked
        self.node_ids.len()+1-self.masked_markers.len()
    }

    fn num_extremities(&self) -> usize {
       if self.degree(TELOMERE).unwrap_or(0) > 0 {
            2*self.num_markers()+1
       } else {
            2*self.num_markers()
       }
    }

    fn trim_singlethread(&mut self, min_size : usize) {
        let nbefore = self.num_markers();
        let to_remove = self.identify_removal_nodes(min_size);
        eprintln!("Identified {} markers for removal.",to_remove.len());
        self.remove_all(&to_remove);
        if self.num_markers() < nbefore/10 {
            eprintln!("Warning: Only {} markers left after trimming",self.num_markers())
        }
    }

    

    fn from_hash_maps(sizes : HashMap<Marker,usize>, adj :HashMap<Extremity,HashSet<Extremity>> ,  nids : HashMap<String,Marker>) -> Self {
        let max_markers = nids.values().max().copied().unwrap_or(0);
        let max_extremities = adj.keys().max().copied().unwrap_or(0);
        let mut node_sizes =vec![0;max_markers+1];
        let mut adjacencies = vec![Vec::new();max_extremities+1];
        for (marker,size) in sizes {
            node_sizes[marker] = size;
        }

        for (x,neighb) in adj {
            adjacencies[x] =  neighb.iter().copied().collect();
            if neighb.contains(&x) {
                adjacencies[x].push(x);
                //selfies.insert(x);
            }
        }

        
        let rnids = reverse_map(&nids);
        let mut masked_markers = HashSet::new();
        //mask markers without id
        for i in 0..(max_markers+1) {
            if !rnids.contains_key(&i) {
                masked_markers.insert(i);
            }
        }


        MBG { node_sizes: node_sizes, adjacencies: adjacencies, node_ids: nids, masked_markers: masked_markers }
    }


    fn from_gfa(path: &str) -> io::Result<Self>{
    if !path.ends_with(".gz") {
        eprintln!("Trying to read uncompressed gfa.");
        let mut rdr = ReaderBuilder::new().has_headers(false).delimiter(b'\t').flexible(true).from_path(path)?;
        return Self::gfa_from_any_reader(&mut rdr);
    } else{
        eprintln!("Trying to read compressed gfa.");
        let fl = File::open(path)?;
        let gz = MultiGzDecoder::new(fl);
        let mut rdr =  ReaderBuilder::new().has_headers(false).delimiter(b'\t').flexible(true).from_reader(gz);
        return Self::gfa_from_any_reader(&mut rdr);
    }
}


fn fill_telomeres(&mut self) {
    let mut new_telos = Vec::new();
    let mmax = self.markers().filter(|x| !self.masked_markers.contains(x)).max().unwrap_or(0);
    let hmax = head(mmax);
    Self::fill_up_vec(&mut self.adjacencies, hmax);
    for (u,adj) in self.adjacencies.iter_mut().enumerate() {
        if u > 1 && !self.masked_markers.contains( &marker(u)) && adj.len() == 0 {
            new_telos.push(u);
            adj.push(TELOMERE);
        }
    }
    self.adjacencies[TELOMERE].extend_from_slice(&new_telos);
}

fn from_unimog(path : &str) -> io::Result<MBG> {
    let node_sizes = Vec::new();
    let mut adjacencies :Vec<Vec<Extremity>> = Vec::new();
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
            
            //eprintln!("{}",m);
            let is_forward;
            let marker;
            (curr_id,is_forward,marker) = parse_marker(&mut node_ids, m, curr_id);
            //eprintln!("is_forward {}",is_forward);
            if first.is_none() {
                first=Some((is_forward,marker));
            }
            if let Some(last) = last {
                let (xta,xtb) = to_adjacency(last,(is_forward,marker));
                Self::fill_up_vec(&mut adjacencies, xta.max(xtb));
                adjacencies[xta].push(xtb);
                adjacencies[xtb].push(xta);
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
            Self::fill_up_vec(&mut adjacencies, xta.max(xtb));
            if line.starts_with(")") {
                adjacencies[xta].push(xtb);
                adjacencies[xtb].push(xta);
            } else if line.starts_with("|") {
                adjacencies[TELOMERE].push(xtb);
                adjacencies[TELOMERE].push(xta);
                adjacencies[xta].push(TELOMERE);
                adjacencies[xtb].push(TELOMERE);
            } else {
                return Err(io::Error::new(io::ErrorKind::Other,"Invalid chromosome end."));
            }
        }

    }
    Ok(MBG { node_sizes: node_sizes, adjacencies: adjacencies, node_ids: node_ids, masked_markers: HashSet::from([TELOMERE])})
}

fn name_to_marker(&self,name : &str) -> Option<Marker> {
    self.node_ids.get(name).filter(|x| !self.masked_markers.contains(*x)).copied()
}

fn marker_names(&self) -> HashMap<Marker,String> {
    reverse_map(&self.node_ids)
}

fn trim_multithread(&mut self, min_size : usize, n_threads : usize) {

    if n_threads <= 1 {
        self.trim_singlethread(min_size);
        return;
    }
    let max_xt = self.adjacencies.len();
    let mut results: Vec<(HashSet<usize>, Vec<(usize,Vec<usize>)>)> = Vec::new();
    let nbefore = self.num_markers();
    thread::scope(|scope| {
        let mut handles = Vec::new();
        let slice_size = max_xt/n_threads +1;
        for i in 0..n_threads {
            let g = &self;
            let lb = (slice_size*i).min(max_xt);
            let rb = (slice_size*(i+1)).min(max_xt);
            eprintln!("Spawning thread {i} processing extremities with index {lb} to {rb} (total {max_xt})");
            let x =  scope.spawn(move || __trim_vertices(g,lb, rb, min_size));
            handles.push(x);
        }
        for x in handles {
            results.push(x.join().unwrap());
        }
        
    });
    eprintln!("Joining results.");
    for (masked, adj) in results {
        self.masked_markers.extend(&masked);
        for m in &masked {
            self.adjacencies[head(*m)] = Vec::new();
            self.adjacencies[tail(*m)] = Vec::new();
        }
        for (i, adjs) in adj {
            self.adjacencies[i]=adjs
        }
    }
    if self.num_markers() < nbefore/10 {
            eprintln!("Warning: Only {} markers left after trimming",self.num_markers())
    }
}
    
    
    
}






pub type Marker = usize;
pub type Extremity = usize;
pub type Adjacency = (Extremity,Extremity);


#[inline(always)]
pub fn head(n : Marker) -> Extremity {
    2*n+1
}

#[inline(always)]
pub fn tail(n : Marker) -> Extremity {
    2*n
}

#[inline(always)]
pub fn other(n : Extremity) -> Extremity{
    if n == TELOMERE {
        return TELOMERE;
    }
    if is_tail(n) {
        head(marker(n))
    } else {
        tail(marker(n))
    }
}

#[inline(always)]
pub fn is_tail(n:Extremity) -> bool{
    match n%2 {
        0 => true,
        _ => false
    }
}

#[inline(always)]
pub fn hdtl_fmt(n :Extremity) -> String {
    let x = match n%2 {
        0 => "t",
        _ => "h"
    };
    marker(n).to_string()+x
}

#[inline(always)]
pub fn marker(xtr : Extremity) -> Marker {
    xtr/2
}

#[cfg(test)]
mod tests {
    use std::fs::read_dir;

    use super::*;
    #[test]
    fn test_hdtl() {
        assert_eq!(4,tail(2));
        assert_eq!(5,head(2));
    }

    #[test]
    fn test_tail() {
        for i in 0..42 {
            assert!(is_tail(tail(i)));
            assert!(!is_tail(head(i)));
        }
    }
    #[test]
    fn test_marker_conversion() {
        for i in 0..42 {
            assert_eq!(i,marker(head(i)));
            assert_eq!(i,marker(tail(i)));
        }
    }

    #[test]
    fn test_trim_graph() {
        mtest_trim_graph::<UBG>(1);
        for i in 1..11 {
            mtest_trim_graph::<MBG>(i);
        }
    }

    fn mtest_trim_graph<T>(n_threads : usize)
    where T : RearrangementGraph
     {
        let mut adj = HashMap :: new();
        let mut node_siz = HashMap :: new();
        let mut node_ids = HashMap :: new();
        node_siz.insert(1, 3);
        node_siz.insert(2, 2);
        node_siz.insert(3, 1);
        for i in 1..4 {
            node_ids.insert(i.to_string(), i);
        }
        for i in 1..=3 {
            adj.insert(head(i), HashSet::new());
            adj.insert(tail(i), HashSet::new());
        }
        adj.get_mut(&head(1)).expect("!").insert(tail(2));
        adj.get_mut(&head(2)).expect("!").insert(tail(3));
        adj.get_mut(&head(3)).expect("!").insert(tail(1));

        adj.get_mut(&tail(1)).expect("!").insert(head(3));
        adj.get_mut(&tail(2)).expect("!").insert(head(1));
        adj.get_mut(&tail(3)).expect("!").insert(head(2));
        let mut ubg = T::from_hash_maps(node_siz, adj, node_ids);


        eprintln!("mtest nthreads {n_threads}");
        if n_threads == 1 {
            ubg.trim_singlethread(2);
        } else {
            ubg.trim_multithread(2, n_threads);
        }
        
        general_ubg_sanity_check(&ubg);
        assert!(ubg.node_size(1).is_some());
        assert!(ubg.node_size(2).is_some());
        assert!(ubg.node_size(3).is_none());

        assert!(ubg.adj_neighbors(head(3)).is_none());
        assert!(ubg.adj_neighbors(tail(3)).is_none());
        //eprintln!("{:?}",ubg.adjacencies);
        let nh2 : HashSet<Extremity> = ubg.adj_neighbors(head(2)).unwrap().collect();
        assert!(nh2.contains(&tail(1)));
        let nt1 : HashSet<Extremity> = ubg.adj_neighbors(tail(1)).unwrap().collect();
        assert!(nt1.contains(&head(2)));
        let mut remaining = HashSet::new();
        for i in 1..=2 {
            remaining.insert(head(i));
            remaining.insert(tail(i));
        }
        for (x,y) in ubg.iter_adjacencies() {
            eprintln!("{x} {y}");
            assert!(remaining.contains(&x));
            assert!(remaining.contains(&y));
        };
    }

    fn general_ubg_sanity_check(ubg: &impl RearrangementGraph) {
        //forbidden telomere thing
        assert!(ubg.degree(1).is_none());
        if ubg.degree(TELOMERE).unwrap_or(0) > 0 {
            assert_eq!(ubg.num_extremities(),ubg.num_markers()*2+1);
        } else {
            assert_eq!(ubg.num_extremities(),ubg.num_markers()*2)
        }
        
        for m in ubg.markers() {
            assert!(ubg.degree(tail(m)).is_some());
            assert!(ubg.degree(head(m)).is_some());
        }
        let mset : HashSet<Marker> = ubg.markers().collect();
        for x in ubg.extremities() {
            if !mset.contains(&marker(x)) && x!=TELOMERE {
                panic!("Extremity {x} exists, but its marker does not!");
            }
            if let Some(neighbors) = ubg.adj_neighbors(x) {
                for y in neighbors {
                    if !mset.contains(&marker(y)) && y!=TELOMERE {
                        panic!("Extremity {y} exists, but its marker does not!");
                    }
                let mut has_x = false;
                for z in ubg.adj_neighbors(y).unwrap() {
                    has_x = has_x || (z==x);
                }
                assert!(has_x);
            }
            }
            
        }
    }

    #[test]
    fn test_read_unimog() {
        mtest_read_unimog::<UBG>();
        mtest_read_unimog::<MBG>();
    }
    fn mtest_read_unimog<T>()
    where T : RearrangementGraph
     {
        let ubg =  T::from_unimog("testfiles/test01.ug").expect("File should be readable");
        general_ubg_sanity_check(&ubg);
        for i in 1..=3 {
            assert!(ubg.adj_neighbors(tail(i)).is_some());
            assert!(ubg.adj_neighbors(head(i)).is_some());
        }
        assert!(ubg.adj_neighbors(head(2)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(1),tail(2)])));
        assert!(ubg.adj_neighbors(tail(2)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(3),head(2)])));
        assert!(ubg.adj_neighbors(head(1)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(2)])));
        assert!(ubg.adj_neighbors(head(3)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(1)])));
        assert!(ubg.adj_neighbors(tail(3)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(2)])));
        assert!(ubg.adj_neighbors(tail(1)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(3)])));
    }

    fn collect_neighbors(ubg : &impl RearrangementGraph,xtr : Extremity) -> HashSet<Extremity> {
        ubg.adj_neighbors(xtr).unwrap().collect()
    }

    #[test]
    fn test_unimog_new() {
        mtest_unimog_new::<UBG>();
        mtest_unimog_new::<MBG>();
    }
    fn mtest_unimog_new<T>()
    where T : RearrangementGraph
     {
        let ubg = T::from_unimog("testfiles/test05.ug").expect("File should be readable");
        general_ubg_sanity_check(&ubg);
        for i in ["1","2","5","6","7"] {
            assert!(ubg.name_to_marker(i).is_some());
        }
        for i in ["3","4","+5","+1"] {
            assert!(ubg.name_to_marker(i).is_none());
        }
        let nd1 = ubg.name_to_marker("1").unwrap();
        let nd2 = ubg.name_to_marker("2").unwrap();
        let nd5 = ubg.name_to_marker("5").unwrap();
        let nd6 = ubg.name_to_marker("6").unwrap();
        let nd7 = ubg.name_to_marker("7").unwrap();
        assert!(collect_neighbors(&ubg,tail(nd1)).eq(&HashSet::from([TELOMERE])));
        assert!(collect_neighbors(&ubg,head(nd1)).eq(&HashSet::from([head(nd2),tail(nd6)])));
        assert!(collect_neighbors(&ubg,head(nd2)).eq(&HashSet::from([head(nd1)])));
        assert!(collect_neighbors(&ubg,tail(nd2)).eq(&HashSet::from([tail(nd5)])));
        assert!(collect_neighbors(&ubg,tail(nd5)).eq(&HashSet::from([tail(nd2)])));
        assert!(collect_neighbors(&ubg,head(nd5)).eq(&HashSet::from([TELOMERE])));
        assert!(collect_neighbors(&ubg,tail(nd7)).eq(&HashSet::from([head(nd7)])));
        assert!(collect_neighbors(&ubg,tail(nd6)).eq(&HashSet::from([head(nd1)])));
        assert!(collect_neighbors(&ubg,head(nd6)).eq(&HashSet::from([TELOMERE])));
    }


    #[test]
    fn test_self_cyle() {
        eprintln!("self cycle test ubg");
        mtest_self_cycle::<UBG>();
        eprintln!("self cycle test mbg");
        mtest_self_cycle::<MBG>();
    }

    fn mtest_self_cycle<T>()
    where T : RearrangementGraph
    {
        let mut sizes = HashMap::new();
        let mut adj = HashMap::new();
        let mut nids = HashMap::new();
        nids.insert(String::from("A"), 1);
        sizes.insert(1,1);
        adj.insert(head(1),HashSet::from([head(1)]));
        adj.insert(tail(1), HashSet::from([TELOMERE]));
        adj.insert(TELOMERE, HashSet::from([tail(1)]));
        let g = T::from_hash_maps(sizes, adj, nids);
        let (a,b) = calc_carp_measure_multithread(&g, 1);
        let mut aexp = Vec::new();
        aexp.push((head(1),head(1)));
        assert_eq!(a,aexp);
    }

    #[test]
    fn test_read_gfa_path() {
        let mut mbg = MBG::from_gfa("testfiles/test08.gfa").unwrap();
        general_ubg_sanity_check(&mbg);
        let tels : HashSet<Extremity> = mbg.adj_neighbors(TELOMERE).unwrap().collect();
        let mut expect : HashSet<Extremity> = HashSet::new();
        let m1 = mbg.name_to_marker(&"1").unwrap();
        let m2 = mbg.name_to_marker(&"2").unwrap();
        expect.insert(tail(m1));
        expect.insert(tail(m2));
        assert_eq!(tels,expect);
        mbg.trim_multithread(2, 3);
        general_ubg_sanity_check(&mbg);
        let m1nb : HashSet<Extremity> = mbg.adj_neighbors(head(m1)).unwrap().collect();
        let mut expectnb = HashSet::new();
        expectnb.insert(TELOMERE);
        expectnb.insert(head(m1));
        assert_eq!(m1nb,expectnb);
        mbg.trim_multithread(4, 10);
        assert_eq!(mbg.num_extremities(),0);
        assert_eq!(mbg.adj_neighbors(TELOMERE).map(|x| x.collect()).unwrap_or(HashSet::new()),HashSet::new());
    }


    #[test]
    fn test_read_gfa_walk() {
        let mut mbg = MBG::from_gfa("testfiles/test10.gfa").unwrap();
        general_ubg_sanity_check(&mbg);
        let tels : HashSet<Extremity> = mbg.adj_neighbors(TELOMERE).unwrap().collect();
        let mut expect : HashSet<Extremity> = HashSet::new();
        let m1 = mbg.name_to_marker(&"1").unwrap();
        let m2 = mbg.name_to_marker(&"2").unwrap();
        expect.insert(tail(m1));
        expect.insert(tail(m2));
        assert_eq!(tels,expect);
        mbg.trim_multithread(2, 3);
        general_ubg_sanity_check(&mbg);
        let m1nb : HashSet<Extremity> = mbg.adj_neighbors(head(m1)).unwrap().collect();
        let mut expectnb = HashSet::new();
        expectnb.insert(TELOMERE);
        expectnb.insert(head(m1));
        assert_eq!(m1nb,expectnb);
        mbg.trim_multithread(4, 10);
        assert_eq!(mbg.num_extremities(),0);
        assert_eq!(mbg.adj_neighbors(TELOMERE).map(|x| x.collect()).unwrap_or(HashSet::new()),HashSet::new());
    }


    #[test]
    fn test_read_gfa_path_single() {
        let mbg = MBG::from_gfa("testfiles/test09.gfa").unwrap();
        general_ubg_sanity_check(&mbg);
        let tels : HashSet<Extremity> = mbg.adj_neighbors(TELOMERE).unwrap().collect();
        let mut expect : HashSet<Extremity> = HashSet::new();
        let m1 = mbg.name_to_marker(&"1").unwrap();
        expect.insert(head(m1));
        expect.insert(tail(m1));
        assert_eq!(tels,expect);
    }


    #[test]
    fn test_read_gfa_walk_single() {
        let mbg = MBG::from_gfa("testfiles/test11.gfa").unwrap();
        general_ubg_sanity_check(&mbg);
        let tels : HashSet<Extremity> = mbg.adj_neighbors(TELOMERE).unwrap().collect();
        let mut expect : HashSet<Extremity> = HashSet::new();
        let m1 = mbg.name_to_marker(&"1").unwrap();
        expect.insert(head(m1));
        expect.insert(tail(m1));
        assert_eq!(tels,expect);
    }

    #[test]
    fn test_random_gfa() {
        for i in read_dir("testfiles/random/").expect("W") {
            let x = i.expect("AAAAAAA");
            let tmpval = x.path();
            let gfafile = tmpval.to_str().unwrap();
            eprintln!("Tesing {gfafile}");
            let mut ubg = UBG::from_gfa(gfafile).expect("X");
            ubg.fill_telomeres();
            let mut mbg = MBG::from_gfa(gfafile).expect("Y");
            mbg.fill_telomeres();
            eprintln!("Sanity checking UBG.");
            general_ubg_sanity_check(&ubg);
            eprintln!("Sanity checking MBG.");
            general_ubg_sanity_check(&mbg);
            eprintln!("Eq check");
            equivalence_check(&ubg, &mbg);
            for flt in 0..10 {
                eprintln!("Filtersize or smth: {flt}");
                let mut ubg_ = ubg.clone();
                let mut mbg_ = mbg.clone();
                ubg_.trim_singlethread(flt);
                mbg_.trim_singlethread(flt);
                
                eprintln!("Sanity checking UBG.");
                general_ubg_sanity_check(&ubg_);
                eprintln!("Sanity checking MBG.");
                general_ubg_sanity_check(&mbg_);
                eprintln!("Eq check");
                equivalence_check(&ubg_, &mbg_);
                for i in 3..5 {
                    let mut mbg_mtt = mbg.clone();
                    mbg_mtt.trim_multithread(flt, i);
                    general_ubg_sanity_check(&mbg_mtt);
                    equivalence_check(&mbg_, &mbg_mtt);
                }
                let (c,uc) = calc_carp_measure_naive(&ubg_);
                let (cc,ucc) = calc_carp_measure_naive(&mbg_);
                let ac = carp_measure_from_adjacencies(&ubg_.iter_adjacencies().collect());
                let acc = carp_measure_from_adjacencies(&mbg_.iter_adjacencies().collect());
                eprintln!("diff ja lol ey {:?}",c.symmetric_difference(&cc));
                assert_eq!(c,cc);
                assert_eq!(uc,ucc);
                assert_eq!(ac,acc);
                assert_eq!(c.len(),ac);
                assert_eq!(cc.len(),acc);
                carp_sanity_check(&mbg_, &c, &uc);
            }

           
        }
    }

    #[test]
    fn test_read_gfa_fail() {
        mtest_read_gfa_fail::<UBG>();
        mtest_read_gfa_fail::<MBG>();
    }
    fn mtest_read_gfa_fail<T>()
    where T : RearrangementGraph
     {
        for path in ["test01.ug","testfiles/test04.gfa","testfiles/test06.gfa","testfiles/test07.gfa"] {
            let x = T::from_gfa(path);
            assert!(x.is_err());
        } 
    }

    #[test]
    fn test_emoji() {
        mtest_emoji::<UBG>();
        mtest_emoji::<MBG>();
    }
    fn mtest_emoji<T>()
    where T : RearrangementGraph
    {
        let mut ubg = T::from_gfa("testfiles/test03.gfa").expect("File should be readable");
        //eprintln!("nodes {:?}",ubg.node_ids);
        //eprintln!("sizes {:?}",ubg.node_sizes);
        ubg.fill_telomeres();
        general_ubg_sanity_check(&ubg);
        //assert!(ubg.adjacencies.get(&tail(1)).expect(".").eq(&HashSet::from([tail(2)])));
        let h1n : HashSet<Extremity> = ubg.adj_neighbors(tail(1)).unwrap().collect();
        assert!(h1n.eq(&HashSet::from([tail(2)])));
        let t2n : HashSet<Extremity> = ubg.adj_neighbors(tail(2)).unwrap().collect();
        assert!(t2n.eq(&HashSet::from([tail(1)])));
        let ttn : HashSet<Extremity> = ubg.adj_neighbors(TELOMERE).unwrap().collect();
        assert!(ttn.eq(&HashSet::from([head(1),head(2)])));
    }

    #[test]
    fn test_read_gfa() {
        mtest_read_gfa::<UBG>();
        mtest_read_gfa::<MBG>();
    }
    fn mtest_read_gfa<T>()
    where T : RearrangementGraph
     {
        let mut ubg = T::from_gfa("testfiles/test02.gfa").expect("File should be readable");
        ubg.fill_telomeres();
        general_ubg_sanity_check(&ubg);
        for i in 1..=4 {
            assert!(ubg.adj_neighbors(tail(i)).is_some());
            assert!(ubg.adj_neighbors(head(i)).is_some());
        }
        let correct_sizes = HashMap::from([(1,5),(2,2),(3,0),(4,6)]);
        for i in 1..=4 {
            assert_eq!(ubg.node_size(i).unwrap(),*correct_sizes.get(&i).expect("."));
        }
        assert!(ubg.adj_neighbors(head(1)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(4)])));
        assert!(ubg.adj_neighbors(head(2)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(4)])));
        assert!(ubg.adj_neighbors(tail(4)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(1),head(2)])));
        assert!(ubg.adj_neighbors(tail(3)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(2)])));
        assert!(ubg.adj_neighbors(tail(2)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(3)])));
        assert!(ubg.adj_neighbors(head(3)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(4)])));
        assert!(ubg.adj_neighbors(tail(1)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(4)])));
        assert!(ubg.adj_neighbors(head(4)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(1),head(3)])));
    }

    fn carp_multithread_sanity_check(ubg : &impl RearrangementGraph, contested : &HashSet<Adjacency>, uncontested : &HashSet<Adjacency>) {
        for i in 1..8 {
            let (c, uc) = calc_carp_measure_multithread(ubg, i);
            assert_eq!(find_dups(&c),Vec::new());
            assert_eq!(find_dups(&uc),Vec::new());
            eprintln!("El difference: {:?}",contested.symmetric_difference(&c.iter().copied().collect()));
            assert_eq!(c.len(),contested.len());
            assert_eq!(uc.len(),uncontested.len());
            let cs : HashSet<Adjacency> = c.iter().copied().collect();
            let ucs : HashSet<Adjacency> = uc.iter().copied().collect();
            assert_eq!(cs,contested.clone());
            assert_eq!(ucs,uncontested.clone());
        }
    }

    fn carp_sanity_check(ubg : &impl RearrangementGraph, contested : &HashSet<Adjacency>, uncontested : &HashSet<Adjacency>) {
        let mut all = HashSet::new();
        for (x,y) in ubg.iter_adjacencies() {
            all.insert(canonicize((x,y)));
            
        }
        //eprintln!("Mief: {:?}",contested.difference(&all).sorted());
        assert!(contested.is_subset(&all));
        assert!(uncontested.is_subset(&all));
        assert!(contested.is_disjoint(uncontested));
        let union : HashSet<Adjacency> = contested.union(uncontested).cloned().collect();
        //eprintln!("Sym diff: {:?}",all.symmetric_difference(&union));
        assert_eq!(all,union);
        for (x,y) in contested {
            assert!(ubg.adj_neighbors(*x).unwrap().collect::<Vec<_>>().len()>1 
                || ubg.adj_neighbors(*y).unwrap().collect::<Vec<_>>().len()>1 || *x==*y);
            assert!(*x!=TELOMERE && *y!=TELOMERE);
        }
        for (x,y) in uncontested {
            assert!((ubg.adj_neighbors(*x).unwrap().collect::<Vec<_>>().len()==1 
                && ubg.adj_neighbors(*y).unwrap().collect::<Vec<_>>().len()==1 && *x!=*y)
                || *x==TELOMERE || *y==TELOMERE);
            assert!(ubg.adj_neighbors(*x).unwrap().collect::<HashSet<_>>().contains(y)
                && ubg.adj_neighbors(*y).unwrap().collect::<HashSet<_>>().contains(x))
        }
        carp_multithread_sanity_check(ubg, contested, uncontested);

    }

    #[test]
    fn test_carp_eva() {
        let ubg = UBG::from_gfa("testfiles/test02.gfa").expect("File should be readable");
        let mbg = MBG::from_gfa("testfiles/test02.gfa").expect("File should be readable");
        let (contested,uncontested) = calc_carp_measure_naive(&ubg);
        let (mcontested,muncontested) = calc_carp_measure_naive(&mbg);
        assert_eq!(contested,mcontested);
        assert_eq!(uncontested,muncontested);
        carp_sanity_check(&ubg, &contested, &uncontested);
        carp_sanity_check(&mbg, &mcontested, &muncontested);
        assert!(mcontested.len()==4);
        assert!(muncontested.len()==1);
    }

    #[test]
    fn test_back_map() {
        let m = HashMap::from([(1,2),(3,4),(5,6)]);
        let r = reverse_map(&m);
        for (a,b) in m.iter() {
            assert!(r.get(b).expect(".")==a)
        }
        for (a,b) in r.iter() {
            assert!(m.get(b).expect(".")==a)
        }
    }


    fn equivalence_check(g1 : &impl RearrangementGraph, g2 : &impl RearrangementGraph) {
        let a1 : HashSet<Adjacency> = g1.iter_adjacencies().map(|x| canonicize(x)).collect();
        let a2 : HashSet<Adjacency> = g2.iter_adjacencies().map(|x| canonicize(x)).collect();
        let mut mp = g1.marker_names();
        let mp2 = g2.marker_names();
        mp.extend(mp2);
        let a_u1 : HashSet<Adjacency> = a1.difference(&a2).copied().collect();
        eprintln!("{a_u1:?}");
        let a_u2 : HashSet<Adjacency> = a2.difference(&a1).copied().collect();
        eprintln!("{a_u2:?}");
        let a_u : HashSet<String> = a1.symmetric_difference(&a2).copied().map(|a| pretty_adjacency(&mp, a)).collect();
        eprintln!("{a_u:?}");
        assert!(a_u.len()==0);
        let m1 : HashSet<Marker> = g1.markers().collect();
        let m2 : HashSet<Marker> = g2.markers().collect();
        eprintln!("{m1:?}");
        assert_eq!(m1,m2);
        let ms1 : HashSet<(Marker,usize)> = m1.iter().map(|x| 
            {
            //eprintln!("{x}");
            (*x,g1.node_size(*x).unwrap())
            }).collect();
        let ms2 : HashSet<(Marker,usize)> = m2.iter().map(|x| (*x,g2.node_size(*x).unwrap())).collect();
        assert_eq!(ms1,ms2);
    }

    #[test]
    fn test_trim_ypestis() {
        let mut mbg = MBG::from_gfa("testfiles/test_ypestis.gfa").expect("File should be readable");
        let mut ubg = UBG::from_gfa("testfiles/test_ypestis.gfa").expect("File should be readable");
        mbg.fill_telomeres();
        ubg.fill_telomeres();
        equivalence_check(&mbg, &ubg);
        eprintln!("Pre-Testing MBG");
        general_ubg_sanity_check(&mbg);
        eprintln!("Pre-Testing UBG");
        general_ubg_sanity_check(&ubg);
        for i in 32..35 {
            let mut ubg = ubg. clone();
            ubg.trim_singlethread(i);
            for nt in 1..9 {
                let mut mbg = mbg. clone();
                mbg.trim_multithread(i,nt);
                eprintln!("Testing MBG t:{i}");
                general_ubg_sanity_check(&mbg);
                eprintln!("Testing UBG t:{i}");
                general_ubg_sanity_check(&ubg);
                equivalence_check(&mbg, &ubg);
            }
            
        }
    }

    #[test]
    fn test_multithread_trimming_empty() {
        let mut mbg = MBG::from_gfa("testfiles/test_16s.gfa").expect("File should be readable");
        mbg.fill_telomeres();
        mbg.trim_multithread(ONE_MILLION,2);
        general_ubg_sanity_check(&mbg);
    }





    #[test]
    fn test_scan() {
        let mut mbg = MBG::from_gfa("testfiles/test12.gfa").unwrap();
        mbg.fill_telomeres();
        general_ubg_sanity_check(&mbg);
        let m1 = mbg.name_to_marker("1").unwrap();
        let m2 = mbg.name_to_marker("2").unwrap();
        let m3 = mbg.name_to_marker("3").unwrap();
        let m4 = mbg.name_to_marker("4").unwrap();
        let m5 = mbg.name_to_marker("5").unwrap();
        let m6 = mbg.name_to_marker("6").unwrap();
        let m7 = mbg.name_to_marker("7").unwrap();
        let m8 = mbg.name_to_marker("8").unwrap();
        let m9 = mbg.name_to_marker("9").unwrap();
        let m10 = mbg.name_to_marker("10").unwrap();

        //Test for emptiness
        let adj = adjacency_neighborhood(m2,24,&mbg);
        assert_eq!(adj,HashSet::new());


        //Test emptiness for self loops
        let adj = adjacency_neighborhood(m1,24,&mbg);
        assert_eq!(adj,HashSet::new());

        //Test complete component
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m2,500,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m2),tail(m4)),
            (head(m4),head(m2)),
            (head(m1),tail(m2)),
            (tail(m1),tail(m1)),
            (head(m3),tail(m1)),
            (tail(m3),tail(m3)),
            (tail(m3),tail(m2)),
            
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test cutoff
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m2,100,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m2),tail(m4)),
            (head(m4),head(m2)),
            (head(m1),tail(m2)),
            (tail(m3),tail(m2))
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test no parallel
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m5,1500,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m5),tail(m6)),
            (head(m5),head(m10)),
            (tail(m5),tail(m10)),
            (head(m9),tail(m10)),
            (tail(m9),tail(m8)),
            (head(m6),TELOMERE),
            (head(m8),TELOMERE)
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test directed in
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m9,1500,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m5),head(m10)),
            (tail(m5),tail(m10)),
            (head(m9),tail(m10)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE)
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test shortened
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m8,200,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m7),tail(m8)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE)
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test shortened
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m8,400,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m7),tail(m8)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE),
            (tail(m7),tail(m6)),
            (head(m9),tail(m10))
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);

        //Test shortened
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m8,650,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m7),tail(m8)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE),
            (tail(m7),tail(m6)),
            (head(m9),tail(m10)),
            (head(m5),head(m10)),
            (head(m6),TELOMERE)
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);

        //Test shortened
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m8,691,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m7),tail(m8)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE),
            (tail(m7),tail(m6)),
            (head(m9),tail(m10)),
            (head(m5),head(m10)),
            (head(m6),TELOMERE),
            (tail(5),tail(10))
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);
    }
}


fn remove_node_from_adj(adjacencies : &mut HashMap<Extremity,HashSet<Extremity>>, xtr : Extremity) {
    let neighbors=adjacencies.get(&xtr).cloned().expect("aa");
    for x in neighbors.iter() {
        let v: &mut HashSet<Extremity> =  adjacencies.get_mut(&x).expect("A");
        v.remove(&xtr);
    }
    adjacencies.remove(&xtr);
}






#[derive(Copy, Clone, Eq, PartialEq)]
struct State {
    cost: usize,
    position: Extremity,
}

// The priority queue depends on `Ord`.
// Explicitly implement the trait so the queue becomes a min-heap
// instead of a max-heap.
impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        // Notice that we flip the ordering on costs.
        // In case of a tie we compare positions - this step is necessary
        // to make implementations of `PartialEq` and `Ord` consistent.
        other.cost.cmp(&self.cost)
            .then_with(|| self.position.cmp(&other.position))
    }
}

// `PartialOrd` needs to be implemented as well.
impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub fn adjacency_neighborhood(m : Marker,max_depth : usize, graph : &impl RearrangementGraph) -> HashSet<Adjacency>{
    let mut visited = BinaryHeap::new();
    let init_dist = graph.node_size(m).unwrap_or(0)/2;
    let mut adjacencies : HashSet<Adjacency> = HashSet::new();
    let mut min_dist = HashMap::new();
    if init_dist > max_depth {
        return  adjacencies;
    }
    visited.push(State {cost : init_dist,position : head(m)});
    visited.push(State {cost : init_dist,position : tail(m)});
    while let Some(State {cost, position}) = visited.pop(){
        if let Some(neighbors) = graph.adj_neighbors(position) {
            for neigbor in neighbors {
                let adj = if position < neigbor {
                        (position,neigbor)  
                    } else {
                        (neigbor,position)
                    };
                adjacencies.insert(adj);
                if neigbor == TELOMERE {
                    //skip telomeres, they're not real connections
                    continue;
                }
                let oend = other(neigbor);
                
                
                let ndist: usize = cost + graph.node_size(marker(neigbor)).unwrap_or(1);//.node_sizes.get(&marker(*neigbor)).unwrap_or(&1);
                if ndist <= max_depth && *min_dist.get(&oend).unwrap_or(&(ndist+1)) > ndist {
                    visited.push(State{cost: ndist,position: oend});
                    min_dist.insert(oend, ndist);
                }
            }
        }
    }
    adjacencies
}


pub fn partial2gfa(ubg : &impl RearrangementGraph, adjacencies : &HashSet<Adjacency>) {
    let nids = ubg.marker_names();
    let mut nodes = HashMap::new();
    let mut linkstrs = Vec::new();
    for (x,y) in adjacencies {
        let m1n = marker(*x);
        let m2n = marker(*y);
        if let Some(m1) = nids.get(&m1n) {
            if let Some(m2) = nids.get(&m2n) {
                nodes.insert(m1n,m1);
                nodes.insert(m2n,m2);
                let orient1 = if is_tail(*x) {
                    "-"
                } else {
                    "+"
                };
                let orient2 = if is_tail(*y) {
                    "+"
                } else {
                    "-"
                };

                linkstrs.push(format!("L\t{m1}\t{orient1}\t{m2}\t{orient2}\t0M"));
            }
        }
        
    }
    for (k,x) in nodes {
        let mut x :String = format!("S\t{x}\t*");
        if let Some(ns) = ubg.node_size(k) {
            x+=&format!("\tLN:i:{ns}");
        } 
        println!("{}",x);
    }
    for lstr in linkstrs {
        println!("{}",lstr)
    }
}

pub fn carp_measure_from_adjacencies(adjacencies : &HashSet<Adjacency>) -> usize {
    let mut carp_measure = 0;
    let mut degrees : HashMap<Extremity, u32> = HashMap::new();
    for (x,y) in adjacencies {
        let xv= degrees.entry(*x).or_insert(0);
        *xv+=1;
        let yv = degrees.entry(*y).or_insert(0);
        *yv+=1;
    }
    for (x,y) in adjacencies {
        if (*degrees.get(x).unwrap() > 1 || *degrees.get(y).unwrap() > 1) && !(*x== TELOMERE) && !(*y == TELOMERE) {
            carp_measure+=1;
        }
    }
    carp_measure
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



fn parse_marker(node_ids: &mut HashMap<String, Marker>, markerstr: &str, curr_id : Marker) -> (Marker,bool,Marker) {
    let mut workslice = markerstr;
    let mut is_forward = true;
    let mut curr_id = curr_id;
    if workslice.starts_with("+") {
        workslice=&workslice[1..];
    }
    if workslice.starts_with("-") {
        workslice=&workslice[1..];
        is_forward=false;
    }
    let m;
    if !node_ids.contains_key(workslice) {
        node_ids.insert(workslice.to_string(),curr_id);
        m=curr_id;
        curr_id+=1;
    } else {
        m=*node_ids.get(workslice).expect("Contains key, but will not give it.");
    }

    (curr_id,is_forward,m)
}


pub fn to_adjacency((ifa,ma):(bool,Marker),(ifb,mb):(bool,Marker)) -> Adjacency{
    //eprintln!("ma {ma} mb {mb}");
    let xta = if ifa {
        head(ma)
    } else {
        tail(ma)
    };
    let xtb = if ifb {
        tail(mb)
    } else {
        head(mb)
    };
    //eprintln!("Adj {} -> {}",hdtl_fmt(xta),hdtl_fmt(xtb));
    (xta,xtb)
}




fn get_or_set_node_id(node_ids: &mut HashMap<String, Marker>, curr_id: Marker, seg_name: String) -> (Marker,Marker){
    let n_id = *node_ids.entry(seg_name).or_insert(curr_id);
    let new_id = if curr_id==n_id{
        curr_id+1
    } else {
        curr_id
    };
    (new_id,n_id)
}

fn canonicize((a,b):Adjacency) -> Adjacency {
    if a < b {
        (a,b)
    } else {
        (b,a)
    }
}

#[inline(always)]
fn naive_hash(x : usize) -> usize{
    (7727*x+7001)%7919
}
#[inline(always)]
fn is_my_adjacency((x,y) : Adjacency) -> bool {
    let hx = naive_hash(x);
    let hy = naive_hash(y);
    hx < hy || (hx == hy && x <= y)
}


fn find_dups<T>(v: &Vec<T>) -> Vec<T>
where 
    T : Eq,
    T : Hash,
    T: Copy
{
    let mut seen = HashSet::new();
    let mut dups = Vec::new();
    for x in v {
        if seen.contains(x) {
            dups.push(*x)
        }
        seen.insert(*x);
    }
    dups
}

pub fn calc_partial_measure(graph : &impl RearrangementGraph, extremities : &[Extremity],_threadnum: usize) -> (Vec<Adjacency>,Vec<Adjacency>) {
    let mut contested = Vec::new();
    let mut uncontested = Vec::new();
    for x in extremities {
        let degx = graph.degree(*x).unwrap();
        if let Some(neighbors) = graph.adj_neighbors(*x) {
            let mut self_seen = false;
            for y in  neighbors{
                let degy = graph.degree(y).unwrap();
                if !is_my_adjacency((*x,y)) {
                    continue;
                }

                if *x==y {
                    if !self_seen {
                        contested.push(canonicize((*x,y)));
                    }
                    self_seen=true;
                } else if *x==TELOMERE || y==TELOMERE {
                    uncontested.push(canonicize((*x,y)));
                } else if degx > 1 || degy > 1 {
                    contested.push(canonicize((*x,y)));
                } else {
                    uncontested.push(canonicize((*x,y)));
                }
            }
        }
        
    }
    assert_eq!(find_dups(&contested),Vec::new());
    assert_eq!(find_dups(&uncontested),Vec::new());
    (contested,uncontested)
}

pub fn calc_carp_measure_multithread(graph : &impl RearrangementGraph, n_threads : usize) -> (Vec<Adjacency>,Vec<Adjacency>) {
    let mut contested = Vec::new();
    let mut uncontested = Vec::new();
    let extremities : Vec<Extremity> = graph.extremities().collect();
    let tot_len = extremities.len();
    let slice_len = (tot_len/n_threads)+1;
    thread::scope(|scope| {
        let mut handles = Vec::new();
        for i in 0..n_threads {
            let lb = (i*slice_len).min(tot_len);
            let rb = ((i+1)*slice_len).min(tot_len);
            let elist = &extremities[lb..rb];
            let handle = scope.spawn(move || {
                calc_partial_measure(graph, elist, i)
            }
            );
            handles.push(handle);
        }
        let mut i=0;
        for handle in handles {
            eprintln!("Joining {i}.");
            let (c,uc) =  handle.join().unwrap();
            i+=1;
            contested.extend(c);
            uncontested.extend(uc);
        }
     });
    (contested,uncontested)
}



pub fn calc_carp_measure_naive(graph : &impl RearrangementGraph) -> (HashSet<Adjacency>,HashSet<Adjacency>){
    let mut contested = HashSet::new();
    let mut uncontested = HashSet::new();
    let adjset : HashSet<Adjacency>= graph.iter_adjacencies().collect();
    for xtr in graph.extremities() {
        if xtr <= 1 {
            for nxtr in graph.adj_neighbors(xtr).unwrap() {
                uncontested.insert((xtr, nxtr));
            }
            continue
        }
        let is_contested=graph.degree(xtr).unwrap()>1;
        //If vertex has degree > 0, it must have neighbors
        let neighb = graph.adj_neighbors(xtr).unwrap();
        for nxtr in neighb {
            if nxtr <= 1 {
                continue;
            }
            let e =canonicize((nxtr,xtr));
            eprintln!("{e:?}");
            assert!(adjset.contains(&e));
            if is_contested {
                contested.insert(e);
                uncontested.remove(&e);
            } else if !contested.contains(&e){
                uncontested.insert(e);
            }
        }
    }
    (contested,uncontested)
}




pub fn reverse_map<K,V>(m : &HashMap<K,V>) -> HashMap<V,K>
where
K : std::hash::Hash+Clone,
V : Eq+std::hash::Hash+Clone
{
    let mut reversed = HashMap::new();
    for (a,b) in m {
        reversed.insert(b.clone(),a.clone());
    }
    reversed
}


fn pretty_extremity(m : &HashMap<Marker,String>,x : Extremity) -> String {
    if !m.contains_key(&marker(x)) {
        panic!("Extremity {x} not found!");
    }
    let s = m.get(&marker(x)).unwrap().clone();
    let ht = if is_tail(x) {
        "t"
    } else {
        "h"
    };
    s+"_"+ht
}

fn pretty_adjacency(m : &HashMap<Marker,String>,(x,y) : Adjacency) -> String{
    pretty_extremity(m, x)+"-"+&pretty_extremity(m, y)
}