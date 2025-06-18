use std::collections::{HashMap};
use std::io::{self};
use clap::{arg, value_parser, ArgGroup, Command};
use scj_carp_rust::*;


fn scan_graph(graph : &UBG,max_depth :usize) -> HashMap<Marker, usize>{
    eprintln!("Scanning graph...");
    let mut node_complexities = HashMap::new();
    let tot_size = graph.node_ids.len();
    let mut i = 1;
    for (m,_) in graph.node_sizes.iter() {
        if *m == 0 {
            continue;
        }
        //eprintln!("Processing node {}",m);
        let adjacencies = adjacency_neighborhood(*m,max_depth, graph);
        let ci = carp_measure_from_adjacencies(&adjacencies);
        node_complexities.insert(*m,ci);
        if i%(tot_size/100) == 0 {
            let percentage = i*100/tot_size;
            eprintln!("Processed {i}/{tot_size} nodes ({percentage}%).");
        }
        i+=1;
    }
    node_complexities
}

fn histogram(node_complexities : &HashMap<Marker,usize>) -> HashMap<usize,usize> {
    let mut hist : HashMap<usize,usize> = HashMap::new();
    for (_, ci) in node_complexities {
        let count = hist.entry(*ci).or_insert(0);
        *count+=1;
    }
    hist
}

fn top_percentile(node_complexities : &HashMap<Marker,usize>,percentile : f64) -> Vec<Marker> {
    let num_nodes = node_complexities.len();
    let hist = histogram(&node_complexities);
    let mut hist_entries : Vec<usize>= hist.keys().cloned().collect();
    hist_entries.sort();
    let mut count = 0;
    let mut thresh_val = 0;
    for e in hist_entries {
        count+=hist.get(&e).unwrap();
        if count as f64 > percentile * num_nodes as f64 {
            thresh_val=e+1;
            break;
        }
    }
    let mut complex_nodes = Vec::new();
    for (u,ci) in node_complexities {
        if *ci > thresh_val {
            complex_nodes.push(*u);
        }
    }
    complex_nodes
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
        .arg(arg!(-c --"context-len" <c>).value_parser(value_parser!(usize)).default_value("10000"))
        .get_matches();
    
    let thresh = *matches.get_one(&"size-thresh").expect("CLI Parsing gone wrong");
    let contextlen = *matches.get_one(&"context-len").expect("CLI Parsing gone wrong");

    eprintln!("Reading graph...");
    let maybe_graph = match (matches.get_one::<String>("gfa")
            , matches.get_one::<String>("unimog")) {
        (Some(gfaf),_) => parse_gfa(gfaf),
        (_,Some(unimog)) =>  parse_unimog(&unimog),
        (_,_) => Err(io::Error::new(io::ErrorKind::Other,"No file specified."))
    };
    let mut graph = maybe_graph.expect("Something went wrong parsing input files");
    eprintln!("Adding telomeres to complete graph.");
    add_telomeres(&mut graph);
    eprintln!("Trimming graph.");
    trim_graph(&mut graph, thresh);
    add_telomeres(&mut graph);
    let node_c =  scan_graph(&graph, contextlen);
    let inv_map = reverse_map(&graph.node_ids);
    for x in top_percentile(&node_c, 0.99999) {
        println!("{}",inv_map.get(&x).unwrap());
    }
}