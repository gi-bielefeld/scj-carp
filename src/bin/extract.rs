use clap::{arg, value_parser, ArgGroup, Command};
use std::io;

use scj_carp_rust::mbg::MBG;
use scj_carp_rust::rearrangement::RearrangementGraph;
use scj_carp_rust::scan::adjacency_neighborhood;
use scj_carp_rust::gfa::partial2gfa;

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
        (Some(gfaf),_) => MBG::from_gfa(gfaf),
        (_,Some(unimog)) =>  MBG::from_unimog(&unimog),
        (_,_) => Err(io::Error::new(io::ErrorKind::Other,"No file specified."))
    };
    let thresh = *matches.get_one(&"size-thresh").expect("CLI Parsing gone wrong");

    let mut graph = maybe_graph.expect("Something went wrong parsing input files");
    eprintln!("Adding telomeres to complete graph.");
    graph.fill_telomeres();
    if thresh > 0 {
        eprintln!("Trimming graph.");
        graph.trim_singlethread(thresh);
    }
    graph.fill_telomeres();
    let start_node : &String = matches.get_one(&"start-node").expect("CLI Parsing gone wrong");
    let max_dist : usize = *matches.get_one(&"max-dist").expect("CLI Parsing gone wrong");
    let marker = graph.name_to_marker(&start_node).expect("Given node is not part of the (trimmed) graph. Make sure that this node id exists and try a lower size threshold.");
    let adjacencies = adjacency_neighborhood(marker, max_dist, &graph);
    partial2gfa(&graph, &adjacencies);
}