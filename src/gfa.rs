use std::collections::{HashMap,HashSet};
use crate::rearrangement::*;
use crate::rearrangement::Marker;

pub const LEN_PREFIX  : &str = "LN:i:";

pub fn get_or_set_node_id(node_ids: &mut HashMap<String, Marker>, curr_id: Marker, seg_name: String) -> (Marker,Marker){
    let n_id = *node_ids.entry(seg_name).or_insert(curr_id);
    let new_id = if curr_id==n_id{
        curr_id+1
    } else {
        curr_id
    };
    (new_id,n_id)
}


pub fn parse_marker(node_ids: &mut HashMap<String, Marker>, markerstr: &str, curr_id : Marker) -> (Marker,bool,Marker) {
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
