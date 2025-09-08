use std::collections::{HashMap};
use std::fs::File;
use std::io::{self, Write};
use std::process::exit;
use clap::{arg, value_parser, ArgGroup, Command};
use scj_carp_rust::rearrangement::*;
use scj_carp_rust::util::*;
use scj_carp_rust::mbg::MBG;
use scj_carp_rust::scan::*;

fn to_gfa_annotated(graph : &impl RearrangementGraph, annotations : &HashMap<Marker,String>, file : &str) -> std::io::Result<()>
{
    let mut file = File::create(file)?;
    let mnames = graph.marker_names();
    for (mid,mname) in &mnames {
        file.write(format!("S\t{mname}\t*").as_bytes())?;
        if let Some(msiz) = graph.node_size(*mid) {
            file.write(format!("\tLN:i:{msiz}").as_bytes())?;
        }
        if let Some(t) = annotations.get(mid) {
            file.write(format!("\t{t}").as_bytes())?;
        }
        file.write("\n".as_bytes())?;
    }
    for (x,y) in graph.iter_adjacencies() {
        let m1n = marker(x);
        let m2n = marker(y);
        if let Some(m1) = mnames.get(&m1n) {
            if let Some(m2) = mnames.get(&m2n) {
                let orient1 = if is_tail(x) {
                    "-"
                } else {
                    "+"
                };
                let orient2 = if is_tail(y) {
                    "+"
                } else {
                    "-"
                };

                file.write(format!("L\t{m1}\t{orient1}\t{m2}\t{orient2}\t0M\n").as_bytes())?;
            }
        }
    }
    Ok(())
}

fn to_heat_html( x :f64) -> String {
    let redval = (255.0*x) as u8;
    let blueval = 110 as u8;
    format!("#{:02x}00{:02x}",redval,blueval)
}

fn write_hist(hist:&HashMap<usize,usize>,path : &str) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for (a,b) in hist {
        write!(&mut file, "{a}\t{b}\n")?;
    }
    Ok(())
}



fn main() {
    //TODO: make struct
    let cmd = Command::new("scj-carp")
        .arg(arg!(-s --"size-thresh" <st> "Size threshold for nodes (nodes of lower sizes are discarded)")
        .value_parser(value_parser!(usize))
        .default_value("0"))
        .arg(arg!(-g --"gfa" <f> "Specify input as GFA file."))
        .arg(arg!(-u --"unimog" <f> "Specify input as unimog file."))
        .group(ArgGroup::new("infile").args(["gfa","unimog"])
                    .required(true))
        .arg(arg!(-c --"context-len" <c>).value_parser(value_parser!(usize)).default_value("500"))
        .arg(arg!(--"colored-gfa" <f> "Output annotated gfa with complexities."))
        .arg(arg!(--"output-histogram" <f> "Output a histogram of complexities."))
        .arg(arg!(--"lower-percentile" <lo> "Output nodes that lie between the lower and higher percentile to standard output.").value_parser(value_parser!(f64)))
        .arg(arg!(--"higher-percentile" <hi> "Output nodes that lie between the lower and higher percentile to standard output.").value_parser(value_parser!(f64)).default_value("1.00"))
        .arg(arg!(-t --"num-threads" <t> "Number of threads to use in the scanning phase. Default: 1.").value_parser(value_parser!(usize)).default_value("1"))
        .arg(arg!(--"ignore-gfa-overlap").num_args(0));
    
    let matches = cmd.get_matches();
    let mut thresh = *matches.get_one(&"size-thresh").expect("CLI Parsing gone wrong");
    let contextlen = *matches.get_one(&"context-len").expect("CLI Parsing gone wrong");
    let n_threads = *matches.get_one(&"num-threads").expect("CLI parsing gone wrong");
    let ignore_gfa_overlap = matches.get_flag(&"ignore-gfa-overlap");
    let is_gfa = matches.get_one::<String>("gfa").is_some();
    let is_unimog = matches.get_one::<String>("unimog").is_some();
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
    eprintln!("{}",CARP_LOGO);
    eprintln!("Reading graph...");
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
        graph.trim_any(thresh,n_threads);
        graph.fill_telomeres();
    }
    let node_c  =  scan_graph_enum_multithread(&graph, contextlen,n_threads);
    let mn = *node_c.values().min().unwrap();
    let mut mx = * node_c.values().max().unwrap();
    if mx == 0 {
        mx =  1;
    }
    let mut colors = HashMap::new();
    for (marker,carpi) in &node_c {
        colors.insert(*marker, format!("CL:z:{}\tcrp:i:{}",to_heat_html(((*carpi-mn) as f64).log2()/(mx as f64).log2()),carpi));
    }
    //find median node
    //node_c.sort_by(|a,b| a.1.cmp(&b.1));
    //let (node,complexity) = node_c.last().unwrap();
    //let backmap  = graph.marker_names();
    //println!("{} {}",backmap.get(&node).unwrap(),complexity)
    if let Some(colorgfapath) = matches.get_one::<String>("colored-gfa") {
        to_gfa_annotated(&graph, &colors, &colorgfapath).expect("Could not write colored gfa file");
    }
    if let Some(histogrampath) = matches.get_one::<String>("output-histogram") {
        let hist = histogram(&node_c);
        write_hist(&hist,histogrampath).expect("Could not write histogram file.");
    }
    if let (Some(lo),Some(hi)) = (matches.get_one::<f64>("lower-percentile"),matches.get_one::<f64>("higher-percentile")) {
        let mmap = graph.marker_names();
        println!("#Node\tSCJ-CARP-measure in env");
        for marker in top_percentile(&node_c, *lo, *hi) {
            let complexity = node_c.get(&marker).unwrap();
            let markerstring = mmap.get(&marker).unwrap();
            println!("{markerstring}\t{complexity}");
        }
    } 
}