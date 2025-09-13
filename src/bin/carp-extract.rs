use clap::{arg, value_parser, ArgGroup, Command};
use std::io;
use std::process::exit;
use scj_carp_rust::mbg::MBG;
use scj_carp_rust::rearrangement::RearrangementGraph;
use scj_carp_rust::scan::adjacency_neighborhood;
use scj_carp_rust::gfa::partial2gfa;
use scj_carp_rust::measure::carp_measure_from_adjacencies;
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
        .arg(arg!(--"ignore-gfa-overlap").num_args(0))
        .group(ArgGroup::new("overlap").args(["ignore-gfa-overlap","size-thresh"]))
        .get_matches();
    let ignore_gfa_overlap = matches.get_flag(&"ignore-gfa-overlap");
    let is_gfa = matches.get_one::<String>("gfa").is_some();
    let is_unimog = matches.get_one::<String>("unimog").is_some();
    let mut thresh = *matches.get_one(&"size-thresh").expect("CLI Parsing gone wrong");

    if !is_gfa && ignore_gfa_overlap {
        eprintln!("Warning: Not a gfa file. Ignoring --ignore-gfa-overlap flag.")
    }
    if is_unimog && thresh > 0 {
        eprintln!("Warning: Unimog files do not support node sizes. Ignoring --size-thresh flag.");
        thresh = 0;
    }
    if thresh > 0 && !ignore_gfa_overlap {
        eprintln!("Error: A gfa graph can only be trimmed with the --ignore-gfa-overlap flag.");
        exit(1);
    }
    let maybe_graph = match (matches.get_one::<String>("gfa")
            , matches.get_one::<String>("unimog")) {
        (Some(gfaf),_) => MBG::from_gfa(gfaf,ignore_gfa_overlap),
        (_,Some(unimog)) =>  MBG::from_unimog(&unimog),
        (_,_) => Err(io::Error::new(io::ErrorKind::Other,"No file specified."))
    };
    
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
    eprintln!("{}",carp_measure_from_adjacencies(&adjacencies));
    partial2gfa(&graph, &adjacencies);
}
