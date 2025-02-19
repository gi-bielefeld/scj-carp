use std::collections::HashSet;
use std::fs::File;
use std::io::{self, Write};
use clap::{arg, value_parser, ArgGroup, Command};
use scj_carp_rust::*;





fn output_ancestral_adj(ubg : &UBG,uncontested: &HashSet<(u32,u32)>,outfile: &mut File) {
    //println!("Writing ancestral adjacencies...");
    let backmap = reverse_map(&ubg.node_ids);
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
        let xm = backmap.get(&marker(*x)).expect("Retranslating went wrong");
        let ym = backmap.get(&marker(*y)).expect("Retranslating went wrong");
        outfile.write(format!("{xm} {xt}\t{ym} {yt}\n").as_bytes()).expect("Could not write ancestral file.");
        
    }
}

fn measure_to_file(p : &str, m : usize) {
    let mut fl = File::create(p).expect("Could not create measure file");
    fl.write_all(format!("Carp index: {}\n",m).as_bytes()).expect("Could not write to measure file");
}

fn main() {
    let matches = Command::new("scj-carp")
        .arg(arg!(-s --"size-thresh" <st> "Size threshold for nodes (nodes of lower sizes are discarded)")
            .value_parser(value_parser!(usize))
            .default_value("0"))
        .arg(arg!(-g --"gfa" <f> "Specify input as GFA file."))
        .arg(arg!(-u --"unimog" <f> "Specify input as unimog file."))
        .group(ArgGroup::new("infile").args(["gfa","unimog"])
                    .required(true))
        .arg(arg!(-a --"write-ancestor" <p> "Path to write ancestral adjacencies to."))
        .arg(arg!(-m --"write-measure" <p> "Path to write the carp measure to."))
        .get_matches();
    
    let thresh = *matches.get_one(&"size-thresh").expect("CLI Parsing gone wrong");

    println!("Reading graph...");
    let maybe_graph = match (matches.get_one::<String>("gfa")
            , matches.get_one::<String>("unimog")) {
        (Some(gfaf),_) => parse_gfa(gfaf),
        (_,Some(unimog)) =>  parse_unimog(&unimog),
        (_,_) => Err(io::Error::new(io::ErrorKind::Other,"No file specified."))
    };
    let mut graph = maybe_graph.expect("Something went wrong parsing input files");
    println!("Adding telomeres to complete graph.");
    add_telomeres(&mut graph);
    println!("Trimming graph.");
    trim_graph(&mut graph, thresh);
    add_telomeres(&mut graph);
    println!("Calculating carp measure.");
    let (contested, uncontested) = calc_carp_measure(&graph);
    let m = contested.len();
    println!("Carp index: {}",m);
    match matches.get_one::<String>("write-measure") {
        None => (),
        Some(p) => measure_to_file(p, m)
    };
    match matches.get_one::<String>("write-ancestor") {
        None => (),
        Some(p) => output_ancestral_adj(&graph, &uncontested,&mut File::create(p).expect("Could not create output file."))
    };
}
