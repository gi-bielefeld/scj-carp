use std::fs::File;
use std::io::{self, Write};
use clap::{arg, value_parser, ArgGroup, Command};
use scj_carp_rust::mbg::*;
use scj_carp_rust::util::*;
use scj_carp_rust::rearrangement::{RearrangementGraph,output_ancestral_adj};
use scj_carp_rust::measure::calc_carp_measure_multithread;

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
    
    let mut thresh = *matches.get_one(&"size-thresh").expect("CLI Parsing gone wrong");
    let threads = *matches.get_one(&"num-threads").expect("CLI Parsing gone wrong");
    let is_unimog = matches.get_one::<String>("unimog").is_some();
    if is_unimog && thresh > 0 {
        eprintln!("Warning: Unimog files do not support node sizes. Ignoring --size-thresh flag.");
        thresh = 0;
    }

    eprintln!("{}",CARP_LOGO);
    println!("Reading graph...");
    let maybe_graph = match (matches.get_one::<String>("gfa")
            , matches.get_one::<String>("unimog")) {
        (Some(gfaf),_) => MBG::from_gfa(gfaf,true),
        (_,Some(unimog)) =>  MBG::from_unimog(&unimog),
        (_,_) => Err(io::Error::new(io::ErrorKind::Other,"No file specified."))
    };
    let mut graph = maybe_graph.expect("Something went wrong parsing input files");
    println!("Adding telomeres to complete graph.");
    graph.fill_telomeres();
    if thresh > 0 {
        println!("Trimming graph.");
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
