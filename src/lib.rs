use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use csv::ReaderBuilder;

#[derive(Debug)]
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

#[cfg(test)]
mod tests {
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
        let mut adj = HashMap :: new();
        let mut node_siz = HashMap :: new();
        node_siz.insert(1, 3);
        node_siz.insert(2, 2);
        node_siz.insert(3, 1);
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
        let mut ubg = UBG {
            adjacencies : adj,
            node_ids : HashMap::new(),
            node_sizes : node_siz

        };
        println!("{:?}",ubg.adjacencies);
        trim_graph(&mut ubg, 2);
        assert!(ubg.node_sizes.contains_key(&1));
        assert!(ubg.node_sizes.contains_key(&2));
        assert!(!ubg.node_sizes.contains_key(&3));
        assert!(!ubg.adjacencies.contains_key(&head(3)));
        assert!(!ubg.adjacencies.contains_key(&tail(3)));
        println!("{:?}",ubg.adjacencies);
        assert!(ubg.adjacencies.get(&head(2)).expect("!").contains(&tail(1)));
        assert!(ubg.adjacencies.get(&tail(1)).expect("!").contains(&head(2)));
        let mut remaining = HashSet::new();
        for i in 1..=2 {
            remaining.insert(head(i));
            remaining.insert(tail(i));
        }
        for (x,neighb) in ubg.adjacencies {
            assert!(remaining.contains(&x));
            for y in neighb {
                assert!(remaining.contains(&y));
            }
        };
    }

    fn general_ubg_sanity_check(ubg: &UBG) {
        //forbidden telomere thing
        assert!(!ubg.adjacencies.contains_key(&1));
        for (m,_) in ubg.node_sizes.iter() {
            assert!(ubg.adjacencies.contains_key(&tail(*m)));
            assert!(ubg.adjacencies.contains_key(&head(*m)));
        }
        for (x,neighb) in ubg.adjacencies.iter() {
            for y in neighb {
                assert!(ubg.adjacencies.get(y).expect(".").contains(x));
            }
        }
    }

    #[test]
    fn test_read_unimog() {
        let ubg = parse_unimog("testfiles/test01.ug").expect("File should be readable");
        general_ubg_sanity_check(&ubg);
        for i in 1..=3 {
            assert!(ubg.adjacencies.contains_key(&tail(i)));
            assert!(ubg.adjacencies.contains_key(&head(i)));
        }
        assert!(ubg.adjacencies.get(&head(2)).expect(".").eq(&HashSet::from([head(1),tail(2)])));
        assert!(ubg.adjacencies.get(&tail(2)).expect(".").eq(&HashSet::from([tail(3),head(2)])));
        assert!(ubg.adjacencies.get(&head(1)).expect(".").eq(&HashSet::from([head(2)])));
        assert!(ubg.adjacencies.get(&head(3)).expect(".").eq(&HashSet::from([tail(1)])));
        assert!(ubg.adjacencies.get(&tail(3)).expect(".").eq(&HashSet::from([tail(2)])));
        assert!(ubg.adjacencies.get(&tail(1)).expect(".").eq(&HashSet::from([head(3)])));
    }

    #[test]
    fn test_unimog_new() {
        let ubg = parse_unimog("testfiles/test05.ug").expect("File should be readable");
        general_ubg_sanity_check(&ubg);
        for i in ["1","2","5","6","7"] {
            assert!(ubg.node_ids.contains_key(i));
        }
        for i in ["3","4","+5","+1"] {
            assert!(!ubg.node_ids.contains_key(i));
        }
        let nd1 = *ubg.node_ids.get("1").expect(".");
        let nd2 = *ubg.node_ids.get("2").expect(".");
        let nd5 = *ubg.node_ids.get("5").expect(".");
        let nd6 = *ubg.node_ids.get("6").expect(".");
        let nd7 = *ubg.node_ids.get("7").expect(".");
        assert!(ubg.adjacencies.get(&tail(nd1)).expect(".").eq(&HashSet::from([0])));
        assert!(ubg.adjacencies.get(&head(nd1)).expect(".").eq(&HashSet::from([head(nd2),tail(nd6)])));
        assert!(ubg.adjacencies.get(&head(nd2)).expect(".").eq(&HashSet::from([head(nd1)])));
        assert!(ubg.adjacencies.get(&tail(nd2)).expect(".").eq(&HashSet::from([tail(nd5)])));
        assert!(ubg.adjacencies.get(&tail(nd5)).expect(".").eq(&HashSet::from([tail(nd2)])));
        assert!(ubg.adjacencies.get(&head(nd5)).expect(".").eq(&HashSet::from([0])));
        assert!(ubg.adjacencies.get(&tail(nd7)).expect(".").eq(&HashSet::from([head(nd7)])));
        assert!(ubg.adjacencies.get(&tail(nd6)).expect(".").eq(&HashSet::from([head(nd1)])));
        assert!(ubg.adjacencies.get(&head(nd6)).expect(".").eq(&HashSet::from([0])));
    }

    #[test]
    fn test_read_gfa_fail() {
        for path in ["test01.ug","testfiles/test04.gfa","testfiles/test06.gfa","testfiles/test07.gfa"] {
            let x = parse_gfa(path);
            println!("{:?}",x);
            assert!(x.is_err());
        } 
    }

    #[test]
    fn test_emoji() {
        let mut ubg = parse_gfa("testfiles/test03.gfa").expect("File should be readable");
        println!("nodes {:?}",ubg.node_ids);
        println!("sizes {:?}",ubg.node_sizes);
        add_telomeres(&mut ubg);
        general_ubg_sanity_check(&ubg);
        println!("{:?}",ubg.adjacencies);
        assert!(ubg.adjacencies.get(&tail(1)).expect(".").eq(&HashSet::from([tail(2)])));
        assert!(ubg.adjacencies.get(&tail(2)).expect(".").eq(&HashSet::from([tail(1)])));
        assert!(ubg.adjacencies.get(&0).expect(".").eq(&HashSet::from([head(1),head(2)])));
    }

    #[test]
    fn test_read_gfa() {
        let mut ubg = parse_gfa("testfiles/test02.gfa").expect("File should be readable");
        add_telomeres(&mut ubg);
        general_ubg_sanity_check(&ubg);
        for i in 1..=4 {
            assert!(ubg.adjacencies.contains_key(&tail(i)));
            assert!(ubg.adjacencies.contains_key(&head(i)));
        }
        let correct_sizes = HashMap::from([(1,5),(2,2),(3,0),(4,6)]);
        for i in 1..=4 {
            assert_eq!(ubg.node_sizes.get(&i).expect("!"),correct_sizes.get(&i).expect("."));
        }
        assert!(ubg.adjacencies.get(&head(1)).expect(".").eq(&HashSet::from([tail(4)])));
        assert!(ubg.adjacencies.get(&head(2)).expect(".").eq(&HashSet::from([tail(4)])));
        assert!(ubg.adjacencies.get(&tail(4)).expect(".").eq(&HashSet::from([head(1),head(2)])));
        assert!(ubg.adjacencies.get(&tail(3)).expect(".").eq(&HashSet::from([tail(2)])));
        assert!(ubg.adjacencies.get(&tail(2)).expect(".").eq(&HashSet::from([tail(3)])));
        assert!(ubg.adjacencies.get(&head(3)).expect(".").eq(&HashSet::from([head(4)])));
        assert!(ubg.adjacencies.get(&tail(1)).expect(".").eq(&HashSet::from([head(4)])));
        assert!(ubg.adjacencies.get(&head(4)).expect(".").eq(&HashSet::from([tail(1),head(3)])));
    }

    fn carp_sanity_check(ubg : &UBG, contested : &HashSet<(u32,u32)>, uncontested : &HashSet<(u32,u32)>) {
        let mut all = HashSet::new();
        for (x,nghb) in ubg.adjacencies.iter() {
            for y in nghb {
                all.insert(canonicize((*x,*y)));
            }
        }
        assert!(contested.is_subset(&all));
        assert!(uncontested.is_subset(&all));
        assert!(contested.is_disjoint(uncontested));
        let union : HashSet<(u32,u32)> = contested.union(uncontested).cloned().collect();
        assert_eq!(all,union);
        for (x,y) in contested {
            assert!(ubg.adjacencies.get(x).expect("!").len()>1 
                || ubg.adjacencies.get(y).expect("!").len()>1);
            assert!(*x!=0 && *y!=0);
        }
        for (x,y) in uncontested {
            assert!(ubg.adjacencies.get(x).expect("!").len()==1 
                && ubg.adjacencies.get(y).expect("!").len()==1);
            assert!(ubg.adjacencies.get(x).expect("!").contains(y)
                && ubg.adjacencies.get(y).expect("!").contains(x))
        }
    }

    #[test]
    fn test_carp_eva() {
        let ubg = parse_gfa("testfiles/test02.gfa").expect("File should be readable");
        let (contested,uncontested) = calc_carp_measure(&ubg);
        carp_sanity_check(&ubg, &contested, &uncontested);
        assert!(contested.len()==4);
        assert!(uncontested.len()==1);
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
}


fn remove_node_from_adj(adjacencies : &mut HashMap<u32,HashSet<u32>>, xtr : u32) {
    let neighbors=adjacencies.get(&xtr).cloned().expect("aa");
    for x in neighbors.iter() {
        let v: &mut HashSet<u32> =  adjacencies.get_mut(&x).expect("A");
        v.remove(&xtr);
    }
    adjacencies.remove(&xtr);
}



pub fn trim_graph(ubg : &mut UBG,threshold : usize) {
    let mut to_remove = Vec::new();
    for (node,sz) in ubg.node_sizes.iter() {
        if *sz < threshold {
            //println!("Removing {} , i.e. {} (hd) {} (tl)",node,head(*node),tail(*node));
            let hd = head(*node);
            let tl = tail(*node);
            let nh = ubg.adjacencies.get(&hd).expect("Assertion violated: Marker extremity not in adjacencies.").clone();
            let nt = ubg.adjacencies.get(&tl).expect("Assertion violated: Marker extremity not in adjacencies.").clone();
            //add adjacencies between the neighboring markers
            for x in nh.iter() {
                for y in nt.iter() {
                    ubg.adjacencies.get_mut(&x).expect("!").insert(*y);
                    ubg.adjacencies.get_mut(&y).expect("!").insert(*x);
                }
            }
            remove_node_from_adj(&mut ubg.adjacencies, hd);
            remove_node_from_adj(&mut ubg.adjacencies, tl);
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

pub fn add_telomeres(ubg : &mut UBG) {
    println!("{:?}",ubg.node_sizes);
    for (node,_) in ubg.node_sizes.clone().iter() {
        println!("{}",node);
        check_add_tel(ubg, head(*node));
        check_add_tel(ubg, tail(*node));

    }
}

pub fn parse_gfa(path: &str) -> io::Result<UBG>{
    println!("Read gfa.");
    let mut node_sizes = HashMap::new();
    let mut adjacencies = HashMap::new(); 
    let mut node_ids: HashMap<String, u32>   = HashMap::new();
    let mut rdr = ReaderBuilder::new().has_headers(false).delimiter(b'\t').flexible(true).from_path(path)?;
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
            let seg_len = match seg_str {
                None => 0,
                Some(y) => y.len()
            };
            let n_id;
            (curr_id,n_id) = get_or_set_node_id(&mut node_ids, curr_id, seg_name);
            if !node_sizes.contains_key(&n_id) {
                node_sizes.insert(n_id, seg_len);
            }
            
            
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
        node_sizes,
        adjacencies,
        node_ids
    })
}


fn parse_marker(node_ids: &mut HashMap<String, u32>, markerstr: &str, curr_id : u32) -> (u32,bool,u32) {
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


//TODO Make constructor etc
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
        let mut first = None;
        let mut last = None;
        while let Some(j) =  line.find(" ") {
            let m = &line[0..j];
            
            println!("{}",m);
            let is_forward;
            let marker;
            (curr_id,is_forward,marker) = parse_marker(&mut node_ids, m, curr_id);
            println!("is_forward {}",is_forward);
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
                adjacencies.entry(0).or_insert(HashSet::new()).insert(xtb);
                adjacencies.entry(0).or_insert(HashSet::new()).insert(xta);
                adjacencies.entry(xta).or_insert(HashSet::new()).insert(0);
                adjacencies.entry(xtb).or_insert(HashSet::new()).insert(0);
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


fn get_or_set_node_id(node_ids: &mut HashMap<String, u32>, curr_id: u32, seg_name: String) -> (u32,u32){
    let n_id = *node_ids.entry(seg_name).or_insert(curr_id);
    let new_id = if curr_id==n_id{
        curr_id+1
    } else {
        curr_id
    };
    (new_id,n_id)
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
    for (xtr,neighb) in &ubg.adjacencies {
        if *xtr <= 1 {
            continue
        }
        let is_contested=neighb.len()>1;
        for nxtr in neighb {
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
    let mut reversed = HashMap::new();
    for (a,b) in m {
        reversed.insert(b.clone(),a.clone());
    }
    reversed
}
