use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use clap::{arg, value_parser, ArgGroup, Command};
use scj_carp_rust::*;





fn output_ancestral_adj(mid2string : &HashMap<Marker,String>,uncontested: &Vec<Adjacency>,outfile: &mut File) {
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

fn measure_to_file(p : &str, m : usize, nmarkers : usize) {
    let mut fl = File::create(p).expect("Could not create measure file");
    fl.write_all(format!("Number of markers: {}\n",nmarkers).as_bytes()).expect("Could not write to measure file");
    fl.write_all(format!("Carp index: {}\n",m).as_bytes()).expect("Could not write to measure file");
}

fn main() {
    //TODO: make struct
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
        .arg(arg!(-t --"num-threads" <t> "Number of threads to use to calculate SCJ CARP index.").value_parser(value_parser!(usize)).default_value("1"))
        .get_matches();
    
    let thresh = *matches.get_one(&"size-thresh").expect("CLI Parsing gone wrong");
    let threads = *matches.get_one(&"num-threads").expect("CLI Parsing gone wrong");
    eprintln!("{}",CARP_LOGO);
    println!("Reading graph...");
    let maybe_graph = match (matches.get_one::<String>("gfa")
            , matches.get_one::<String>("unimog")) {
        (Some(gfaf),_) => MBG::from_gfa(gfaf),
        (_,Some(unimog)) =>  MBG::from_unimog(&unimog),
        (_,_) => Err(io::Error::new(io::ErrorKind::Other,"No file specified."))
    };
    let mut graph = maybe_graph.expect("Something went wrong parsing input files");
    println!("Adding telomeres to complete graph.");
    graph.fill_telomeres();
    println!("Trimming graph.");
    if thresh > 0 {
        graph.trim_multithread(thresh,threads);
        graph.fill_telomeres();
    }
    println!("Calculating carp measure.");
    let (contested, uncontested) = calc_carp_measure_multithread(&graph,threads);
    let m = contested.len();
    println!("Carp index: {}",m);
    print!("On {} markers",graph.num_markers());
    if let Some(p)=  matches.get_one::<String>("write-measure") {
        measure_to_file(p, m,graph.num_markers());
    }
    if let Some(p)=  matches.get_one::<String>("write-ancestor") {
        output_ancestral_adj(&graph.marker_names(), &uncontested,&mut File::create(p).expect("Could not create output file."));
    }
}
