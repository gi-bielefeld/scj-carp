use std::io::{self, BufRead, BufReader};
use csv::{ReaderBuilder,Reader};
use std::collections::{HashMap, HashSet};
use std::thread;
use flate2::read::MultiGzDecoder;
use itertools::Itertools;
use std::fs::File;

use crate::rearrangement::*;
use crate::util::*;
use crate::gfa::*;

#[derive(Debug,Clone)]
pub struct MBG {
    node_sizes : Vec<usize>,
    adjacencies : Vec<Vec<Extremity>>,
    node_ids : HashMap<String,Marker>,
    masked_markers : HashSet<Marker>,
    ext_overlap : usize
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


    fn gfa_from_any_reader<R>(rdr : &mut Reader<R>,ignore_overlap : bool) -> Result<MBG, io::Error>
    where 
        R : std::io::Read
    {
        //eprintln!("Read gfa.");
        let mut node_sizes = Vec::new();
        let mut adjacencies : Vec<Vec<Extremity>> = Vec::new(); 
        let mut node_ids: HashMap<String, Marker>   = HashMap::new();
        let mut curr_id = 1;
        let mut i :usize = 0;
        let mut n_edges :usize = 0;
        let mut telomeres : HashSet<(String,bool)> = HashSet::new();
        let mut seen_edges = HashSet::new();
        let mut x_overlap : Option<usize> = None;
        let mut warned_overlaps = false;
        let mut warned_cigar = None;
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
                    fill_up_vec(&mut node_sizes, n_id);
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
                    fill_up_vec(&mut adjacencies, axtr);
                }
                if adjacencies.len() <= bxtr {
                    fill_up_vec(&mut adjacencies, bxtr);
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
                if !ignore_overlap {
                    if let Some(cigar) = x.get(5) {
                        let prs_ovlp = cigar.strip_suffix("M").ok_or("Cigar string does not end with M.".to_string()).and_then(|cstr| cstr.parse::<usize>().map_err(|e| e.to_string()));
                        match prs_ovlp {
                            Ok(overlap) => if overlap!=x_overlap.unwrap_or(overlap) {
                                if SAFE_GFA_OVERLAP {
                                    panic!("Not supported: Non-fixed overlaps in Cigar-string. Previous overlap {x_overlap:?}. Current overlap: {overlap}");
                                } else {
                                    warned_overlaps=true;
                                    x_overlap=Some(0);
                                }
                            } else {
                                x_overlap = Some(overlap);
                            },
                            Err(errmsg) => if SAFE_GFA_OVERLAP {
                                panic!("Not supported: Non-match Cigar string overlaps: {cigar}. Error: {errmsg}. If you want to ignore this (at your own risk), recompile with SAFE_GFA_OVERLAP=false.");
                            } else {
                                warned_cigar=Some((cigar.to_owned(),errmsg));
                            }
                        }
                    }
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
                    fill_up_vec(&mut adjacencies, xtr);
            }
            adjacencies.get_mut(xtr).unwrap().push(TELOMERE);
            adjacencies.get_mut(TELOMERE).unwrap().push(xtr);
        }

        let cigproblem = warned_cigar.is_some();
        if let Some((cigar,errmsg)) = warned_cigar {
            eprintln!("Warning: Unsupported Cigar string overlap: {cigar}. Error: {errmsg}.");
        }

        if warned_overlaps {
            eprintln!("Warning: Unsupported overlaps in Cigar strings: Non-fixed length.")
        }

        if warned_overlaps || cigproblem {
            eprintln!("Warning: all overlaps have been set to 0.")
        }
        
        Ok(MBG { node_sizes: node_sizes, adjacencies: adjacencies, node_ids: node_ids, masked_markers: HashSet::from([TELOMERE]) , ext_overlap : x_overlap.unwrap_or(0)})
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
            let nmarkers =  self.num_markers()+1;
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
            if self.ext_overlap > 0 {
                if  SAFE_GFA_OVERLAP {
                    panic!("Cannot safely trim graph with overlaps.");
                } else {
                    eprintln!("Trimming requires 0-Overlap between segments. Overlap will be ignored from now on.");
                    self.ext_overlap=0;
                }
            }
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
        if self.ext_overlap > 0 {
            if  SAFE_GFA_OVERLAP {
                panic!("Cannot safely trim graph with overlaps.");
            } else {
                eprintln!("Trimming requires 0-Overlap between segments. Overlap will be ignored from now on.");
                self.ext_overlap=0;
            }
        }
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


        MBG { node_sizes: node_sizes, adjacencies: adjacencies, node_ids: nids, masked_markers: masked_markers, ext_overlap :0 }
    }


    fn from_gfa(path: &str, ignore_overlap : bool) -> io::Result<Self>{
    if !path.ends_with(".gz") {
        eprintln!("Trying to read uncompressed gfa.");
        let mut rdr = ReaderBuilder::new().has_headers(false).delimiter(b'\t').flexible(true).from_path(path)?;
        return Self::gfa_from_any_reader(&mut rdr,ignore_overlap);
    } else{
        eprintln!("Trying to read compressed gfa.");
        let fl = File::open(path)?;
        let gz = MultiGzDecoder::new(fl);
        let mut rdr =  ReaderBuilder::new().has_headers(false).delimiter(b'\t').flexible(true).from_reader(gz);
        return Self::gfa_from_any_reader(&mut rdr,ignore_overlap);
    }
}


fn fill_telomeres(&mut self) {
    let mut new_telos = Vec::new();
    let mmax = self.markers().filter(|x| !self.masked_markers.contains(x)).max().unwrap_or(0);
    let hmax = head(mmax);
    fill_up_vec(&mut self.adjacencies, hmax);
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
    let mut seen_edges = HashSet::new();
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
                let canone = canonicize((xta,xtb));
                insert_adj(&mut adjacencies, &mut seen_edges, canone);   
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
            fill_up_vec(&mut adjacencies, xta.max(xtb));
            if line.starts_with(")") {
                //adjacencies[xta].push(xtb);
                //adjacencies[xtb].push(xta);
                insert_adj(&mut adjacencies, &mut seen_edges, (xta,xtb));
            } else if line.starts_with("|") {
                //adjacencies[TELOMERE].push(xtb);
                //adjacencies[TELOMERE].push(xta);
                //adjacencies[xta].push(TELOMERE);
                //adjacencies[xtb].push(TELOMERE);
                insert_adj(&mut adjacencies, &mut seen_edges, (TELOMERE,xta));
                insert_adj(&mut adjacencies, &mut seen_edges, (TELOMERE,xtb));
            } else {
                return Err(io::Error::new(io::ErrorKind::Other,"Invalid chromosome end."));
            }
        }

    }
    //eprintln!("nids : {node_ids:?}");
    Ok(MBG { node_sizes: node_sizes, adjacencies: adjacencies, node_ids: node_ids, masked_markers: HashSet::from([TELOMERE]), ext_overlap:0})
}

fn name_to_marker(&self,name : &str) -> Option<Marker> {
    self.node_ids.get(name).filter(|x| !self.masked_markers.contains(*x)).copied()
}

fn marker_names(&self) -> HashMap<Marker,String> {
    reverse_map(&self.node_ids)
}

fn trim_multithread(&mut self, min_size : usize, n_threads : usize) {
    if self.ext_overlap > 0 {
            if  SAFE_GFA_OVERLAP {
                panic!("Cannot safely trim graph with overlaps.");
            } else {
                eprintln!("Trimming requires 0-Overlap between segments. Overlap will be ignored from now on.");
                self.ext_overlap=0;
            }
    }
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
    fn overlap(&self, _:Extremity,_:Extremity) -> usize {
        return self.ext_overlap;
    }
    
    
}

#[inline(always)]
fn insert_adj(adjacencies: &mut Vec<Vec<usize>>, seen_edges: &mut HashSet<(usize, usize)>, adj: Adjacency) {
    let canonic_adj = canonicize(adj);
    let (xta,xtb) = canonic_adj;
    if !seen_edges.contains(&canonic_adj) {
        fill_up_vec(adjacencies, xta.max(xtb));
        adjacencies[xta].push(xtb);
        adjacencies[xtb].push(xta);
        seen_edges.insert(canonic_adj);
    }
}


fn fill_up_vec<T>(adj : &mut Vec<T>, to : usize)
    where T: std::default::Default
     {
        for _ in adj.len()..(to+1) {
            adj.push(T::default());
        }
    }

