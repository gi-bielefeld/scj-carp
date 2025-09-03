use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::thread;
use crate::rearrangement::*;
use crate::measure::*;
use crate::util::*;

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

/* 
fn scan_partial(graph : &impl RearrangementGraph, max_depth :usize , markers : &[usize],thread_num : usize) -> HashMap<Marker,usize> {
    let mut node_complexities = HashMap::new();
    let tot_size = markers.len();
    let about_ten_percent = (markers.len()/100 +1).max(10000);
    //eprintln!("Thread {thread_num}: {about_ten_percent}");
    let mut i = 0;
    for m in markers {
        if *m == 0 {
            continue;
        }
        //eprintln!("Processing node {}",m);
        let adjacencies = adjacency_neighborhood(*m,max_depth, graph);
        let ci = carp_measure_from_adjacencies(&adjacencies);
        node_complexities.insert(*m,ci);
        i+=1;
        if i%(about_ten_percent) == 0 {
            let percentage = i*100/tot_size;
            eprintln!("Thread {thread_num} processed {i}/{tot_size} nodes ({percentage}%).");
        }
        
    }
    node_complexities
}



fn scan_graph_multithread(graph : &impl RearrangementGraph,max_depth :usize, n_threads : usize) -> HashMap<Marker, usize> {
    let markerlist : Vec<Marker> = graph.markers().collect();
    let mut node_complexities = HashMap::new();
    
    let totlen = markerlist.len();
    thread::scope(|scope| {
        let mut handles = Vec::new();
        let slice_size = markerlist.len()/n_threads +1;
        for i in 0..n_threads {
            let lb = slice_size*i;
            let rb = (slice_size*(i+1)).min(markerlist.len());
            eprintln!("Spawning thread {i} processing markers with index {lb} to {rb} (total {totlen})");
            let mlist = &markerlist[lb..rb];
            let x =  scope.spawn(move || scan_partial(graph, max_depth, mlist, i));
            handles.push(x);
        }
        for x in handles {
            let hm = x.join().unwrap();
            node_complexities.extend(hm);
        }
        
    });
    
    node_complexities
}
*/


fn scan_enumerate(graph : &impl RearrangementGraph, max_depth :usize , start : usize, end : usize,thread_num : usize) -> HashMap<Marker,usize> {
    let mut node_complexities = HashMap::new();
    let tot_size = end - start;
    let about_one_percent = (tot_size/100 +1).min(ONE_MILLION);
    //eprintln!("Thread {thread_num}: {about_ten_percent}");
    let mut i = 0;
    for m in start..end {
        if m == 0 || graph.adj_neighbors(head(m)).is_none(){
            continue;
        }
        //eprintln!("Processing node {}",m);
        let adjacencies = adjacency_neighborhood(m,max_depth, graph);
        let ci = carp_measure_from_adjacencies(&adjacencies);
        node_complexities.insert(m,ci);
        i+=1;
        if i%(about_one_percent) == 0 {
            let percentage = i*100/tot_size;
            eprintln!("Thread {thread_num} processed {i}/{tot_size} nodes ({percentage}%).");
        }
        
    }
    node_complexities
}

pub fn scan_graph_enum_multithread(graph : &impl RearrangementGraph,max_depth :usize, n_threads : usize) -> HashMap<Marker, usize> {
    let mmax : Marker = graph.markers().max().unwrap_or(0)+1;
    let mut node_complexities = HashMap::new();
    let totlen = mmax;
    thread::scope(|scope| {
        let mut handles = Vec::new();
        let slice_size = totlen/n_threads +1;
        for i in 0..n_threads {
            let lb = slice_size*i;
            let rb = (slice_size*(i+1)).min(totlen);
            eprintln!("Spawning thread {i} processing markers with index {lb} to {rb} (total {totlen})");
            let x =  scope.spawn(move || scan_enumerate(graph,max_depth,lb,rb,i));
            handles.push(x);
        }
        for x in handles {
            let hm = x.join().unwrap();
            node_complexities.extend(hm);
        }
        
    });
    
    node_complexities
}


pub fn scan_graph(graph : &impl RearrangementGraph,max_depth :usize) -> HashMap<Marker, usize>{
    eprintln!("Scanning graph...");
    let mut node_complexities = HashMap::new();
    let tot_size = graph.num_markers();
    let mut i = 1;
    let mut about_one_percent = tot_size/100;
    if about_one_percent == 0 {
            about_one_percent=1;
    }
    let markerlist : Vec<Marker> = graph.markers().collect();
    for m in &markerlist {
        if *m == 0 {
            continue;
        }
        //eprintln!("Processing node {}",m);
        let adjacencies = adjacency_neighborhood(*m,max_depth, graph);
        let ci = carp_measure_from_adjacencies(&adjacencies);
        node_complexities.insert(*m,ci);
        
        if i%(about_one_percent) == 0 {
            let percentage = i*100/tot_size;
            eprintln!("Processed {i}/{tot_size} nodes ({percentage}%).");
        }
        i+=1;
    }
    node_complexities
}



pub fn histogram(node_complexities : &HashMap<Marker,usize>) -> HashMap<usize,usize> {
    let mut hist : HashMap<usize,usize> = HashMap::new();
    for (_, ci) in node_complexities {
        let count = hist.entry(*ci).or_insert(0);
        *count+=1;
    }
    hist
}

pub fn top_percentile(node_complexities : &HashMap<Marker,usize>,percentile_low : f64,percentile_high : f64) -> Vec<Marker> {
    let num_nodes = node_complexities.len();
    let hist = histogram(&node_complexities);
    let mut hist_entries : Vec<usize>= hist.keys().cloned().collect();
    hist_entries.sort();
    let mut count = 0;
    let mut thresh_low = None;
    let mut thresh_high = None;
    for e in hist_entries {
        let high_enough = (count as f64) >= percentile_low * num_nodes as f64;
        let not_too_high = (count as f64) < percentile_high * num_nodes as f64;
        count+=hist.get(&e).unwrap();
        if high_enough && not_too_high  {
            if thresh_low.is_none() {
                thresh_low = Some(e);
            }
            thresh_high = Some(e+1);
        } else if !not_too_high  {
            break;
        }
    }
    eprintln!("Lowthresh {thresh_low:?} Highthresh {thresh_high:?}");
    let mut complex_nodes = Vec::new();
    if let (Some(lt),Some(ht)) = (thresh_low,thresh_high) {
        for (u,ci) in node_complexities {
        if *ci < ht && *ci >= lt {
            complex_nodes.push(*u);
        }
    }
    }
    
    
    complex_nodes
}