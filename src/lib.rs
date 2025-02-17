use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use csv::ReaderBuilder;

pub struct UBG {
    pub node_sizes : HashMap<u32,usize>,
    pub adjacencies : HashMap<u32,HashSet<u32>>,
    pub node_ids : HashMap<String,u32>
}

pub fn head(n : u32) -> u32 {
    2*n+1
}

pub fn tail(n : u32) -> u32 {
    2*n
}

pub fn is_tail(n:u32) -> bool{
    match n%2 {
        0 => true,
        _ => false
    }
}

pub fn hdtl_fmt(n :u32) -> String {
    let x = match n%2 {
        0 => "t",
        _ => "h"
    };
    marker(n).to_string()+x
}

pub fn marker(xtr : u32) -> u32 {
    xtr/2
}

fn remove_node_from_adj(adjacencies : &mut HashMap<u32,HashSet<u32>>, xtr : u32) {
    let neighbors=adjacencies.get(&xtr).cloned().expect("aa");
    for x in neighbors.iter() {
        let v: &mut HashSet<u32> =  adjacencies.get_mut(&x).expect("A");
        v.remove(&xtr);
    }
    adjacencies.remove(&xtr);
    for x in neighbors.iter() {
        for y in neighbors.iter() {
            if *x!=*y && *x!=xtr && *y != xtr {
                adjacencies.get_mut(x).expect("Neighbor does not have any adjacencies").insert(*y);
            }
        }
    }
}

pub fn trim_graph(ubg : &mut UBG,threshold : usize) {
    let mut to_remove = Vec::new();
    for (node,sz) in ubg.node_sizes.iter() {
        if *sz < threshold {
            //println!("Removing {} , i.e. {} (hd) {} (tl)",node,head(*node),tail(*node));
            remove_node_from_adj(&mut ubg.adjacencies, head(*node));
            remove_node_from_adj(&mut ubg.adjacencies, tail(*node));
            to_remove.push(*node);   
        }
    }

    for node in to_remove {
        ubg.node_sizes.remove(&node);
    }
}

fn check_add_tel(ubg : &mut UBG, n : u32) {
    match ubg.adjacencies.get(&n) {
        None => ubg.adjacencies.insert(n, HashSet::new()),
        Some(_) => None
    };
    match ubg.adjacencies.get(&0) {
        None => ubg.adjacencies.insert(0, HashSet::new()),
        Some(_) => None
    };
    let neighbors= ubg.adjacencies.get_mut(&n).expect("Just inserted value disappeared");
    if neighbors.is_empty() {
        neighbors.insert(0);
        ubg.adjacencies.get_mut(&0).expect("test").insert(n);
    }
    
}

pub fn add_telomoeres(ubg : &mut UBG) {
    for (node,_) in ubg.node_sizes.clone().iter() {
        check_add_tel(ubg, head(*node));
        check_add_tel(ubg, tail(*node));

    }
}

pub fn parse_gfa(path: &str) -> io::Result<UBG>{
    println!("Read gfa.");
    let mut node_sizes = HashMap::new();
    let mut adjacencies = HashMap::new(); 
    let mut node_ids: HashMap<String, u32>   = HashMap::new();
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').flexible(true).from_path(path)?;
    let mut curr_id = 1;
    let mut i :u32 = 0;
    let mut n_edges :usize = 0;
    for res in rdr.records() {
        let x = res?;
        i+=1;
        if i%1000000==0 {
            println!("Read {} lines.",i);
            println!("{} nodes and {} edges in graph.",node_sizes.len(),n_edges);
        }
        if x.len()==0 {
            continue;
        }
        let entrytype = x.get(0);
        let entrytype= match entrytype {
            None => return Err(io::Error::new(io::ErrorKind::Other,"Non-empty line is empty")),
            Some(y) => y
        };
        if entrytype.eq("S") {
            let seg_name = x.get(1);
            let seg_str = x.get(2);
            let seg_name= match seg_name {
                None => return Err(io::Error::new(io::ErrorKind::Other,"Empty segment label")),
                Some(y) => y.to_string()
            };
            let seg_len = match seg_str {
                None => 0,
                Some(y) => y.len()
            };
            let n_id;
            (curr_id,n_id) = get_or_set_node_id(&mut node_ids, curr_id, seg_name);
            if !node_sizes.contains_key(&n_id) {
                node_sizes.insert(n_id, seg_len);
            }
            
            
        } else if entrytype.eq("L") {
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
    
    
    Ok(UBG {
        node_sizes:node_sizes,
        adjacencies:adjacencies,
        node_ids:node_ids
    })
}


pub fn parse_marker(node_ids: &mut HashMap<String, u32>, markerstr: &str, curr_id : u32) -> (u32,bool,u32) {
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


pub fn to_adjacency((ifa,ma):(bool,u32),(ifb,mb):(bool,u32)) -> (u32,u32){
    println!("ma {ma} mb {mb}");
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
    println!("Adj {} -> {}",hdtl_fmt(xta),hdtl_fmt(xtb));
    (xta,xtb)
}

pub fn parse_unimog(path : &str) -> io::Result<UBG> {
    let node_sizes = HashMap::new();
    let mut adjacencies :HashMap<u32, HashSet<u32>> = HashMap::new();
    let mut node_ids: HashMap<String, u32>   = HashMap::new();
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut curr_id = 1;
    for line in reader.lines() {
        let line  = &line?;
        if line.starts_with(&">") {
            continue;
        }
        let mut line = &line[..];
        let mut i = line.find(" ");
        let mut first = None;
        let mut last = None;
        while !i.is_none() {
            let j = i.expect("Suddenly none");
            let m = &line[0..j];
            
            println!("{}",m);
            let is_forward;
            let marker;
            (curr_id,is_forward,marker) = parse_marker(&mut node_ids, m, curr_id);
            println!("is_forward {}",is_forward);
            if first.is_none() {
                first=Some((is_forward,marker));
            }
            if !last.is_none() {
                let (xta,xtb) = to_adjacency(last.expect("there"),(is_forward,marker));
                if !adjacencies.contains_key(&xta) {
                    adjacencies.insert(xta, HashSet::new());
                }
                if !adjacencies.contains_key(&xtb) {
                    adjacencies.insert(xtb, HashSet::new());
                }
                adjacencies.get_mut(&xta).expect("A").insert(xtb);
                adjacencies.get_mut(&xtb).expect("A").insert(xta);
            }
            line = &line[j+1..];
            last = Some((is_forward,marker));
            if line.starts_with(")") || line.starts_with("|") {
                break;
            }
            i = line.find(" ");
        }
        if let Some((ifa,ma)) = last {
            let (ifb,mb) = first.expect("If a last marker exists, there must be a first marker");
            let (xta,xtb) = to_adjacency((ifa,ma),(ifb,mb));
            if !adjacencies.contains_key(&xta) {
                adjacencies.insert(xta, HashSet::new());
            }
            if !adjacencies.contains_key(&xtb) {
                adjacencies.insert(xtb, HashSet::new());
            }
            if line.starts_with(")") {
                adjacencies.get_mut(&xta).expect("A").insert(xtb);
                adjacencies.get_mut(&xtb).expect("A").insert(xta);
            } else if line.starts_with("|") {
                if !adjacencies.contains_key(&0) {
                    adjacencies.insert(0, HashSet::new());
                }
                adjacencies.get_mut(&0).expect("A").insert(xtb);
                adjacencies.get_mut(&0).expect("A").insert(xta);
                adjacencies.get_mut(&xta).expect("A").insert(0);
                adjacencies.get_mut(&xtb).expect("A").insert(0);
            } else {
                return Err(io::Error::new(io::ErrorKind::Other,"Invalid chromosome end."));
            }
        }

    }
    Ok(UBG {
        node_sizes:node_sizes,
        adjacencies:adjacencies,
        node_ids:node_ids
    })
}


fn get_or_set_node_id(node_ids: &mut HashMap<String, u32>, curr_id: u32, seg_name: String) -> (u32,u32){
    let mut new_id = curr_id;
    let n_id;
    if !(node_ids.contains_key(&seg_name)) {
        node_ids.insert(seg_name, curr_id);
        n_id = new_id;
        new_id+=1;
    } else {
        n_id = *node_ids.get(&seg_name).expect("This shouldn't happen.");
    }
    return (new_id,n_id)
}

fn canonicize((a,b):(u32,u32)) -> (u32,u32) {
    if a < b {
        (a,b)
    } else {
        (b,a)
    }
}


pub fn calc_carp_measure(ubg : &UBG) -> (HashSet<(u32,u32)>,HashSet<(u32,u32)>){
    let mut contested = HashSet::new();
    let mut uncontested = HashSet::new();
    for (xtr,neighb) in ubg.adjacencies.iter() {
        if *xtr <= 1 {
            continue
        }
    
        let is_contested=neighb.len()>1;
        for nxtr in neighb.iter() {
            if *nxtr <= 1 {
                continue;
            }
            let e =canonicize((*nxtr,*xtr));

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
    let mut rm = HashMap::new();
    for (a,b) in m {
        
        rm.insert(b.clone(),a.clone());
    }
    rm
}
