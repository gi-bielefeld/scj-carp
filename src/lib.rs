use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::hash::Hash;
use std::io::{self, BufRead, BufReader};
use csv::ReaderBuilder;
use std::cmp::Ordering;
use std::thread;



const ONE_MILLION :usize = 1000000;
const LEN_PREFIX  : &str = "LN:i:";

pub const CARP_LOGO : &str = 
"'kl               You're running SCJ CARP
XMM.                        Version 0.0.1
.MMW.                                    
 MMd                                     
 XM.         Feel free to report bugs at:
 lK      github.com/gi-bielefeld/scj-carp
 .Mx                                     
 .MMO.                                   
  MMMXc.                                 
  KMMMMM0l.                              
  .WMMMMMMM0l.                           
   :MMMMMMMMMM0c.                        
    ;WMMMMMMMMMMMXd;                     
     'NMMMMMMMMMMMMMWOl.                 
      .XMMMMMMMMMMMMMMMMXc               
       .xMMMMMMMMMMMMMMMMMXc             
         ,NMMMMMMMMMMMMMKkxOXO:          
          'NMMMMMMMMMMX,....'NMK'        
        .KMMMMMMMMMMMMN:...,dMMMW'       
        :MMMMMMMMMMMMMMMNXNMMMMMMW.      
         ,okOOo'.;kWMMMMMMMMMMMMMMN.     
                .ckNWNMMMMMMMMMMMMk,.    
             .dOxc,.  .:dxkO0XNNo.       
            ''.                ;         
";
pub const CARP_VERSION : &str = "0.0.1";

#[derive(Debug)]
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
    fn trim(&mut self, min_size : usize);
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
        Some(self.adjacencies.get(&n)?.len())
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

    fn trim(&mut self, threshold : usize) {
        let mut to_remove = Vec::new();
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
            remove_node_from_adj(&mut self.adjacencies, hd);
            remove_node_from_adj(&mut self.adjacencies, tl);
            to_remove.push(*node);   
        }
    }

    for node in to_remove {
        self.node_sizes.remove(&node);
    }
    }


    fn iter_adjacencies(&self) -> impl Iterator<Item=Adjacency> {
        self.adjacencies.iter().flat_map(|(x,neighbors)| {
            neighbors.iter().filter_map(|y| {
                if *y >*x {
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
            eprintln!("is_forward {}",is_forward);
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

fn name_to_marker(&self,name : &str) -> Option<Marker> {
    self.node_ids.get(name).copied()
}


fn marker_names(&self) -> HashMap<Marker,String> {
    reverse_map(&self.node_ids)
}
}


#[derive(Debug)]
pub struct MBG {
    node_sizes : Vec<usize>,
    adjacencies : Vec<Vec<Extremity>>,
    node_ids : HashMap<String,Marker>,
    masked_markers : HashSet<Marker>
}


impl MBG {
    #[inline(always)]
    fn remove_marker(&mut self, m : Marker) {
        if self.masked_markers.contains(&m) {
            return
        }
        let tailnb  = self.adjacencies.get(tail(m)).map(|x: &Vec<usize>|x.clone()).unwrap();
        let headnb  = self.adjacencies.get(head(m)).map(|x|x.clone()).unwrap();
        self.adjacencies[head(m)]= Vec::new();
        self.adjacencies[tail(m)] = Vec::new();
        //It's a bit annoying going with the HashMap intermediary
        for x in &tailnb {
            let mut xneighbors : HashSet<Extremity> = self.adjacencies.get(*x).unwrap().iter().copied().collect();
            xneighbors.remove(&tail(m));
            for y in &headnb {
                xneighbors.insert(*y);
            }
            self.adjacencies[*x] = xneighbors.iter().copied().collect();
        }

        for x in &headnb {
            let mut xneighbors : HashSet<Extremity> = self.adjacencies.get(*x).unwrap().iter().copied().collect();
            xneighbors.remove(&head(m));
            for y in &tailnb {
                xneighbors.insert(*y);
            }
            self.adjacencies[*x]=xneighbors.iter().copied().collect();
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



}

impl RearrangementGraph for MBG {
    fn degree(&self,n:Extremity) -> Option<usize> {
        if self.masked_markers.contains(&marker(n)) && n!=0 {
            return None;
        }
        self.adjacencies.get(n).map(|x| x.len())
    }

    fn adj_neighbors(&self,n:Extremity) -> Option<impl Iterator<Item=Extremity>> {
        if self.masked_markers.contains(&marker(n)) && n!=0 {
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
        self.markers().flat_map(|m| [tail(m),head(m)]).filter(|x| !self.masked_markers.contains(&marker(*x))).chain([0])
    }

    fn iter_adjacencies(&self) -> impl Iterator<Item=Adjacency> {
        self.adjacencies.iter().enumerate().flat_map(|(xtr,neighb)|{
         neighb.iter().filter_map(move |ytr| {
        if xtr < *ytr {
            Some((xtr,*ytr))
        } else {
            None
            }})    
        })
    }

    fn node_size(&self,n:Marker) -> Option<usize> {
        if n == 0 || self.masked_markers.contains(&n) {
            return None;
        }
        self.node_sizes.get(n).copied()
    }

    fn num_markers(&self) -> usize {
        //+1 because of telomere node, which does not have an id, but occurs in masked
        self.node_ids.len()-self.masked_markers.len()+1
    }

    fn num_extremities(&self) -> usize {
       if self.degree(0).unwrap_or(0) > 0 {
            2*self.num_markers()+1
       } else {
            2*self.num_markers()
       }
    }

    fn trim(&mut self, min_size : usize) {
        let mut to_remove = Vec::new();
        let nbefore = self.num_markers();
        for (m,size) in self.node_sizes.iter().enumerate() {
            if *size < min_size {
                to_remove.push(m);
            }
        }
        for m in to_remove {
            self.remove_marker(m);
        }
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
    eprintln!("Read gfa.");
    let mut node_sizes = Vec::new();
    let mut adjacencies : Vec<Vec<Extremity>> = Vec::new(); 
    let mut node_ids: HashMap<String, Marker>   = HashMap::new();
    let mut rdr = ReaderBuilder::new().has_headers(false).delimiter(b'\t').flexible(true).from_path(path)?;
    let mut curr_id = 1;
    let mut i :usize = 0;
    let mut n_edges :usize = 0;
    
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
            adjacencies.get_mut(axtr).unwrap().push(bxtr);
            adjacencies.get_mut(bxtr).unwrap().push(axtr);
            n_edges+=1;
        }

    }
    //eprintln!("{node_sizes:?}");
    
    Ok(MBG { node_sizes: node_sizes, adjacencies: adjacencies, node_ids: node_ids, masked_markers: HashSet::from([0]) })
}


fn fill_telomeres(&mut self) {
    let mut new_telos = Vec::new();
    let mmax = self.markers().filter(|x| !self.masked_markers.contains(x)).max().unwrap_or(0);
    let hmax = head(mmax);
    Self::fill_up_vec(&mut self.adjacencies, hmax);
    for (u,adj) in self.adjacencies.iter_mut().enumerate() {
        if u > 1 && !self.masked_markers.contains( &marker(u)) && adj.len() == 0 {
            new_telos.push(u);
            adj.push(0);
        }
    }
    self.adjacencies[0].extend_from_slice(&new_telos);
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
            
            eprintln!("{}",m);
            let is_forward;
            let marker;
            (curr_id,is_forward,marker) = parse_marker(&mut node_ids, m, curr_id);
            eprintln!("is_forward {}",is_forward);
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
                adjacencies[0].push(xtb);
                adjacencies[0].push(xta);
                adjacencies[xta].push(0);
                adjacencies[xtb].push(0);
            } else {
                return Err(io::Error::new(io::ErrorKind::Other,"Invalid chromosome end."));
            }
        }

    }
    Ok(MBG { node_sizes: node_sizes, adjacencies: adjacencies, node_ids: node_ids, masked_markers: HashSet::from([0])})
}

fn name_to_marker(&self,name : &str) -> Option<Marker> {
    self.node_ids.get(name).filter(|x| !self.masked_markers.contains(*x)).copied()
}

fn marker_names(&self) -> HashMap<Marker,String> {
    reverse_map(&self.node_ids)
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
    if n == 0 {
        return 0;
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
        mtest_trim_graph::<UBG>();
        mtest_trim_graph::<MBG>();
    }

    fn mtest_trim_graph<T>()
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

        ubg.trim(2);
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
        for m in ubg.markers() {
            assert!(ubg.degree(tail(m)).is_some());
            assert!(ubg.degree(head(m)).is_some());
        }
        for x in ubg.extremities() {
            if let Some(neighbors) = ubg.adj_neighbors(x) {
                for y in neighbors {
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
        assert!(collect_neighbors(&ubg,tail(nd1)).eq(&HashSet::from([0])));
        assert!(collect_neighbors(&ubg,head(nd1)).eq(&HashSet::from([head(nd2),tail(nd6)])));
        assert!(collect_neighbors(&ubg,head(nd2)).eq(&HashSet::from([head(nd1)])));
        assert!(collect_neighbors(&ubg,tail(nd2)).eq(&HashSet::from([tail(nd5)])));
        assert!(collect_neighbors(&ubg,tail(nd5)).eq(&HashSet::from([tail(nd2)])));
        assert!(collect_neighbors(&ubg,head(nd5)).eq(&HashSet::from([0])));
        assert!(collect_neighbors(&ubg,tail(nd7)).eq(&HashSet::from([head(nd7)])));
        assert!(collect_neighbors(&ubg,tail(nd6)).eq(&HashSet::from([head(nd1)])));
        assert!(collect_neighbors(&ubg,head(nd6)).eq(&HashSet::from([0])));
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
        let ttn : HashSet<Extremity> = ubg.adj_neighbors(0).unwrap().collect();
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

    fn carp_sanity_check(ubg : &impl RearrangementGraph, contested : &HashSet<Adjacency>, uncontested : &HashSet<Adjacency>) {
        let mut all = HashSet::new();
        for (x,y) in ubg.iter_adjacencies() {
            all.insert(canonicize((x,y)));
            
        }
        assert!(contested.is_subset(&all));
        assert!(uncontested.is_subset(&all));
        assert!(contested.is_disjoint(uncontested));
        let union : HashSet<Adjacency> = contested.union(uncontested).cloned().collect();
        assert_eq!(all,union);
        for (x,y) in contested {
            assert!(ubg.adj_neighbors(*x).unwrap().collect::<Vec<_>>().len()>1 
                || ubg.adj_neighbors(*y).unwrap().collect::<Vec<_>>().len()>1);
            assert!(*x!=0 && *y!=0);
        }
        for (x,y) in uncontested {
            assert!(ubg.adj_neighbors(*x).unwrap().collect::<Vec<_>>().len()==1 
                && ubg.adj_neighbors(*y).unwrap().collect::<Vec<_>>().len()==1);
            assert!(ubg.adj_neighbors(*x).unwrap().collect::<HashSet<_>>().contains(y)
                && ubg.adj_neighbors(*y).unwrap().collect::<HashSet<_>>().contains(x))
        }
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
                if neigbor == 0 {
                    //skip telomeres, they're not real connections
                    continue;
                }
                let oend = other(neigbor);
                let adj = if position < neigbor {
                  (position,neigbor)  
                } else {
                    (neigbor,position)
                };
                adjacencies.insert(adj);
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
        if *degrees.get(x).unwrap() > 1 || *degrees.get(y).unwrap() > 1 {
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
    eprintln!("ma {ma} mb {mb}");
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
    eprintln!("Adj {} -> {}",hdtl_fmt(xta),hdtl_fmt(xtb));
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

pub fn calc_partial_measure(graph : &impl RearrangementGraph, extremities : &[Extremity],threadnum: usize) -> (Vec<Adjacency>,Vec<Adjacency>) {
    let mut contested = Vec::new();
    let mut uncontested = Vec::new();
    for x in extremities {
        let degx = graph.degree(*x).unwrap();
        if let Some(neighbors) = graph.adj_neighbors(*x) {
            for y in  neighbors{
                let degy = graph.degree(y).unwrap();
                if !is_my_adjacency((*x,y)) {
                    continue;
                }
                if degx > 1 || degy > 1 {
                    contested.push(canonicize((*x,y)));
                } else {
                    uncontested.push(canonicize((*x,y)));
                }
            }
        }
        
    }

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
            let lb = i*slice_len;
            let rb = ((i+1)*slice_len).min(tot_len);
            let elist = &extremities[lb..rb];
            let handle = scope.spawn(move || {
                calc_partial_measure(graph, elist, i)
            }
            );
            handles.push(handle);
        }
        for handle in handles {
            let (c,uc) =  handle.join().unwrap();
            contested.extend(c);
            uncontested.extend(uc);
        }
     });
    (contested,uncontested)
}



pub fn calc_carp_measure_naive(graph : &impl RearrangementGraph) -> (HashSet<Adjacency>,HashSet<Adjacency>){
    let mut contested = HashSet::new();
    let mut uncontested = HashSet::new();
    for xtr in graph.extremities() {
        if xtr <= 1 {
            continue
        }
        let is_contested=graph.degree(xtr).unwrap_or(0)>1;
        //If vertex has degree > 0, it must have neighbors
        let neighb = graph.adj_neighbors(xtr).unwrap();
        for nxtr in neighb {
            if nxtr <= 1 {
                continue;
            }
            let e =canonicize((nxtr,xtr));
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
