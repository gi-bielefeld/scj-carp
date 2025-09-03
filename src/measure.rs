use std::collections::{HashMap, HashSet};
use std::thread;
use crate::rearrangement::*;
use crate::util::*;


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

