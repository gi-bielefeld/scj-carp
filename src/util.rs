use std::collections::{HashMap,HashSet};
use std::hash::Hash;
use crate::rearrangement::*;


pub const ONE_MILLION :usize = 1000000;

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


pub fn pretty_extremity(m : &HashMap<Marker,String>,x : Extremity) -> String {
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

pub fn pretty_adjacency(m : &HashMap<Marker,String>,(x,y) : Adjacency) -> String{
    pretty_extremity(m, x)+"-"+&pretty_extremity(m, y)
}

pub fn find_dups<T>(v: &Vec<T>) -> Vec<T>
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


#[inline(always)]
pub fn naive_hash(x : usize) -> usize{
    (7727*x+7001)%7919
}
#[inline(always)]
pub fn is_my_adjacency((x,y) : Adjacency) -> bool {
    let hx = naive_hash(x);
    let hy = naive_hash(y);
    hx < hy || (hx == hy && x <= y)
}