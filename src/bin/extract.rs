use clap::{arg, value_parser, ArgGroup, Command};
use scj_carp_rust::{parse_gfa,parse_unimog,add_telomeres,adjacency_neighborhood,partial2gfa,trim_graph};
use std::io;

fn main() {
    let matches = Command::new("filter")
        .arg(arg!(-s --"size-thresh" <st> "Size threshold for nodes (nodes of lower sizes are discarded)")
            .value_parser(value_parser!(usize))
            .default_value("0"))
        .arg(arg!(-g --"gfa" <f> "Specify input as GFA file."))
        .arg(arg!(-u --"unimog" <f> "Specify input as unimog file."))
        .group(ArgGroup::new("infile").args(["gfa","unimog"])
                    .required(true))
        .arg(arg!(-n --"start-node"<n> "start extracting from node"))
        .arg(arg!(-d --"max-dist" <d> "Maximum distance from start node").value_parser(value_parser!(usize)).required(true))
        .get_matches();
    let maybe_graph = match (matches.get_one::<String>("gfa")
            , matches.get_one::<String>("unimog")) {
        (Some(gfaf),_) => parse_gfa(gfaf),
        (_,Some(unimog)) =>  parse_unimog(&unimog),
        (_,_) => Err(io::Error::new(io::ErrorKind::Other,"No file specified."))
    };
    let thresh = *matches.get_one(&"size-thresh").expect("CLI Parsing gone wrong");

    let mut graph = maybe_graph.expect("Something went wrong parsing input files");
    eprintln!("Adding telomeres to complete graph.");
    add_telomeres(&mut graph);
    if thresh > 0 {
        eprintln!("Trimming graph.");
        trim_graph(&mut graph, thresh);
    }
    add_telomeres(&mut graph);
    let start_node : &String = matches.get_one(&"start-node").expect("CLI Parsing gone wrong");
    let max_dist : usize = *matches.get_one(&"max-dist").expect("CLI Parsing gone wrong");
    let marker = *graph.node_ids.get(start_node).unwrap();
    let adjacencies = adjacency_neighborhood(marker, max_dist, &graph);
    partial2gfa(&graph, &adjacencies);
}