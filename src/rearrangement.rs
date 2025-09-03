use std::collections::{HashMap, HashSet};
use std::io;
use std::fs::File;
use std::io::Write;

pub type Marker = usize;
pub type Extremity = usize;
pub type Adjacency = (Extremity,Extremity);

pub const TELOMERE : Extremity = 0;


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

pub fn canonicize((a,b):Adjacency) -> Adjacency {
    if a < b {
        (a,b)
    } else {
        (b,a)
    }
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


pub fn output_ancestral_adj(mid2string : &HashMap<Marker,String>,uncontested: &Vec<Adjacency>,outfile: &mut File) {
    //println!("Writing ancestral adjacencies...");
    for (x,y) in uncontested {
        let xt = match is_tail(*x) {
            true => "t",
            false => "h"
        };
        let yt = match is_tail(*y){
            true => "t",
            false => "h"
        };
        if *x == 0 || *y==0 {
            continue;
        }
        let xm = mid2string.get(&marker(*x)).expect("Retranslating went wrong");
        let ym = mid2string.get(&marker(*y)).expect("Retranslating went wrong");
        outfile.write(format!("{xm} {xt}\t{ym} {yt}\n").as_bytes()).expect("Could not write ancestral file.");
        
    }
}