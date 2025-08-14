use std::collections::{HashMap};
use std::fs::File;
use std::io::{self, Write};
use clap::{arg, value_parser, ArgGroup, Command};
use scj_carp_rust::*;
use std::thread;

const ONE_MILLION :usize = 1000000;

fn scan_partial(graph : &impl RearrangementGraph, max_depth :usize , markers : &[usize],thread_num : usize) -> HashMap<Marker,usize> {
    let mut node_complexities = HashMap::new();
    let tot_size = markers.len();
    let about_ten_percent = (markers.len()/100 +1);//.max(10000);
    //eprintln!("Thread {thread_num}: {about_ten_percent}");
    let mut i = 0;
    for m in markers {
        if *m == 0 {
            continue;
        }
        //eprintln!("Processing node {}",m);
        let adjacencies = adjacency_neighborhood(*m,max_depth, graph);
        let ci = carp_measure_from_adjacencies(&adjacencies);
        node_complexities.insert(*m,ci);
        i+=1;
        if i%(about_ten_percent) == 0 {
            let percentage = i*100/tot_size;
            eprintln!("Thread {thread_num} processed {i}/{tot_size} nodes ({percentage}%).");
        }
        
    }
    node_complexities
}



fn scan_enumerate(graph : &impl RearrangementGraph, max_depth :usize , start : usize, end : usize,thread_num : usize) -> HashMap<Marker,usize> {
    let mut node_complexities = HashMap::new();
    let tot_size = end - start;
    let about_five_percent = (tot_size/20 +1).min(ONE_MILLION);
    //eprintln!("Thread {thread_num}: {about_ten_percent}");
    let mut i = 0;
    for m in start..end {
        if m == 0 || graph.adj_neighbors(head(m)).is_none(){
            continue;
        }
        //eprintln!("Processing node {}",m);
        let adjacencies = adjacency_neighborhood(m,max_depth, graph);
        let ci = carp_measure_from_adjacencies(&adjacencies);
        node_complexities.insert(m,ci);
        i+=1;
        if i%(about_five_percent) == 0 {
            let percentage = i*100/tot_size;
            eprintln!("Thread {thread_num} processed {i}/{tot_size} nodes ({percentage}%).");
        }
        
    }
    node_complexities
}



fn scan_graph_multithread(graph : &impl RearrangementGraph,max_depth :usize, n_threads : usize) -> HashMap<Marker, usize> {
    let markerlist : Vec<Marker> = graph.markers().collect();
    let mut node_complexities = HashMap::new();
    
    let totlen = markerlist.len();
    thread::scope(|scope| {
        let mut handles = Vec::new();
        let slice_size = markerlist.len()/n_threads +1;
        for i in 0..n_threads {
            let lb = slice_size*i;
            let rb = (slice_size*(i+1)).min(markerlist.len());
            eprintln!("Spawning thread {i} processing markers with index {lb} to {rb} (total {totlen})");
            let mlist = &markerlist[lb..rb];
            let x =  scope.spawn(move || scan_partial(graph, max_depth, mlist, i));
            handles.push(x);
        }
        for x in handles {
            let hm = x.join().unwrap();
            node_complexities.extend(hm);
        }
        
    });
    
    node_complexities
}



fn scan_graph_enum_multithread(graph : &impl RearrangementGraph,max_depth :usize, n_threads : usize) -> HashMap<Marker, usize> {
    let mmax : Marker = graph.markers().max().unwrap_or(0)+1;
    let mut node_complexities = HashMap::new();
    
    let totlen = mmax - 1;
    thread::scope(|scope| {
        let mut handles = Vec::new();
        let slice_size = totlen/n_threads +1;
        for i in 0..n_threads {
            let lb = slice_size*i;
            let rb = (slice_size*(i+1)).min(totlen);
            eprintln!("Spawning thread {i} processing markers with index {lb} to {rb} (total {totlen})");
            let x =  scope.spawn(move || scan_enumerate(graph,max_depth,lb,rb,i));
            handles.push(x);
        }
        for x in handles {
            let hm = x.join().unwrap();
            node_complexities.extend(hm);
        }
        
    });
    
    node_complexities
}


fn scan_graph(graph : &impl RearrangementGraph,max_depth :usize) -> HashMap<Marker, usize>{
    eprintln!("Scanning graph...");
    let mut node_complexities = HashMap::new();
    let tot_size = graph.num_markers();
    let mut i = 1;
    let mut about_one_percent = tot_size/100;
    if about_one_percent == 0 {
            about_one_percent=1;
    }
    let markerlist : Vec<Marker> = graph.markers().collect();
    for m in &markerlist {
        if *m == 0 {
            continue;
        }
        //eprintln!("Processing node {}",m);
        let adjacencies = adjacency_neighborhood(*m,max_depth, graph);
        let ci = carp_measure_from_adjacencies(&adjacencies);
        node_complexities.insert(*m,ci);
        
        if i%(about_one_percent) == 0 {
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

fn top_percentile(node_complexities : &HashMap<Marker,usize>,percentile_low : f64,percentile_high : f64) -> Vec<Marker> {
    let num_nodes = node_complexities.len();
    let hist = histogram(&node_complexities);
    let mut hist_entries : Vec<usize>= hist.keys().cloned().collect();
    hist_entries.sort();
    let mut count = 0;
    let mut thresh_low = None;
    let mut thresh_high = None;
    for e in hist_entries {
        let high_enough = count as f64 > percentile_low * num_nodes as f64 ;
        let not_too_high = count as f64 <= percentile_high * num_nodes as f64;
        count+=hist.get(&e).unwrap();
        if high_enough && not_too_high  {
            if thresh_low.is_none() {
                thresh_low = Some(e);
            }
            thresh_high = Some(e+1);
        } else if !not_too_high  {
            break;
        }
    }
    eprintln!("{thresh_low:?} {thresh_high:?}");
    let mut complex_nodes = Vec::new();
    if let (Some(lt),Some(ht)) = (thresh_low,thresh_high) {
        for (u,ci) in node_complexities {
        if *ci < ht && *ci >= lt {
            complex_nodes.push(*u);
        }
    }
    }
    
    
    complex_nodes
}

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
    let blueval = (255.0*0.5*(1.0-x)) as u8;
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
        .arg(arg!(-c --"context-len" <c>).value_parser(value_parser!(usize)).default_value("10000"))
        .arg(arg!(--"colored-gfa" <f> "Output annotated gfa with complexities."))
        .arg(arg!(--"output-histogram" <f> "Output a histogram of complexities."))
        .arg(arg!(--"lower-percentile" <lo> "Output nodes that lie between the lower and higher percentile to standard output. Default 0.99").value_parser(value_parser!(f64)).default_value("0.99"))
        .arg(arg!(--"higher-percentile" <hi> "Output nodes that lie between the lower and higher percentile to standard output. Default 1.00").value_parser(value_parser!(f64)).default_value("1.00"))
        .arg(arg!(-t --"num-threads" <t> "Number of threads to use in the scanning phase. Default: 1.").value_parser(value_parser!(usize)).default_value("1"));
    
    let matches = cmd.get_matches();
    let thresh = *matches.get_one(&"size-thresh").expect("CLI Parsing gone wrong");
    let contextlen = *matches.get_one(&"context-len").expect("CLI Parsing gone wrong");
    let n_threads = *matches.get_one(&"num-threads").expect("CLI parsing gone wrong");
    eprintln!("{}",CARP_LOGO);
    eprintln!("Reading graph...");
    let maybe_graph = match (matches.get_one::<String>("gfa")
            , matches.get_one::<String>("unimog")) {
        (Some(gfaf),_) => MBG::from_gfa(gfaf),
        (_,Some(unimog)) =>  MBG::from_unimog(&unimog),
        (_,_) => Err(io::Error::new(io::ErrorKind::Other,"No file specified."))
    };
    let mut graph = maybe_graph.expect("Something went wrong parsing input files");
    eprintln!("Adding telomeres to complete graph.");
    graph.fill_telomeres();
    eprintln!("Trimming graph.");
    graph.trim_multithread(thresh,n_threads);
    graph.fill_telomeres();
    let node_c  =  scan_graph_enum_multithread(&graph, contextlen,n_threads);
    let mn = *node_c.values().min().unwrap();
    let mut mx = * node_c.values().max().unwrap();
    if mx == 0 {
        mx =  1;
    }
    let mut colors = HashMap::new();
    for (marker,carpi) in &node_c {
        colors.insert(*marker, format!("CL:z:{}\tcrp:i:{}",to_heat_html((*carpi-mn) as f64/mx as f64),carpi));
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
        for marker in top_percentile(&node_c, *lo, *hi) {
            let complexity = node_c.get(&marker).unwrap();
            println!("{marker}\t{complexity}");
        }
    } 
}