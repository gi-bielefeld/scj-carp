use std::collections::{HashMap, HashSet};
use std::fs::read_dir;
use crate::util::*;
use crate::rearrangement::*;
use crate::measure::*;
use crate::scan::*;
use crate::ubg::*;
use crate::mbg::*;

    #[test]
    fn test_hdtl() {
        assert_eq!(4,tail(2));
        assert_eq!(5,head(2));
    }

    #[test]
    fn test_tail() {
        for i in 0..42 {
            assert!(is_tail(tail(i)));
            assert!(!is_tail(head(i)));
        }
    }
    #[test]
    fn test_marker_conversion() {
        for i in 0..42 {
            assert_eq!(i,marker(head(i)));
            assert_eq!(i,marker(tail(i)));
        }
    }

    #[test]
    fn test_trim_graph() {
        mtest_trim_graph::<UBG>(1);
        for i in 1..11 {
            mtest_trim_graph::<MBG>(i);
        }
    }

    fn mtest_trim_graph<T>(n_threads : usize)
    where T : RearrangementGraph
     {
        let mut adj = HashMap :: new();
        let mut node_siz = HashMap :: new();
        let mut node_ids = HashMap :: new();
        node_siz.insert(1, 3);
        node_siz.insert(2, 2);
        node_siz.insert(3, 1);
        for i in 1..4 {
            node_ids.insert(i.to_string(), i);
        }
        for i in 1..=3 {
            adj.insert(head(i), HashSet::new());
            adj.insert(tail(i), HashSet::new());
        }
        adj.get_mut(&head(1)).expect("!").insert(tail(2));
        adj.get_mut(&head(2)).expect("!").insert(tail(3));
        adj.get_mut(&head(3)).expect("!").insert(tail(1));

        adj.get_mut(&tail(1)).expect("!").insert(head(3));
        adj.get_mut(&tail(2)).expect("!").insert(head(1));
        adj.get_mut(&tail(3)).expect("!").insert(head(2));
        let mut ubg = T::from_hash_maps(node_siz, adj, node_ids);


        eprintln!("mtest nthreads {n_threads}");
        if n_threads == 1 {
            ubg.trim_singlethread(2);
        } else {
            ubg.trim_multithread(2, n_threads);
        }
        
        general_ubg_sanity_check(&ubg);
        assert!(ubg.node_size(1).is_some());
        assert!(ubg.node_size(2).is_some());
        assert!(ubg.node_size(3).is_none());

        assert!(ubg.adj_neighbors(head(3)).is_none());
        assert!(ubg.adj_neighbors(tail(3)).is_none());
        //eprintln!("{:?}",ubg.adjacencies);
        let nh2 : HashSet<Extremity> = ubg.adj_neighbors(head(2)).unwrap().collect();
        assert!(nh2.contains(&tail(1)));
        let nt1 : HashSet<Extremity> = ubg.adj_neighbors(tail(1)).unwrap().collect();
        assert!(nt1.contains(&head(2)));
        let mut remaining = HashSet::new();
        for i in 1..=2 {
            remaining.insert(head(i));
            remaining.insert(tail(i));
        }
        for (x,y) in ubg.iter_adjacencies() {
            eprintln!("{x} {y}");
            assert!(remaining.contains(&x));
            assert!(remaining.contains(&y));
        };
    }

    fn general_ubg_sanity_check(ubg: &impl RearrangementGraph) {
        //forbidden telomere thing
        assert!(ubg.degree(1).is_none());
        if ubg.degree(TELOMERE).unwrap_or(0) > 0 {
            assert_eq!(ubg.num_extremities(),ubg.num_markers()*2+1);
        } else {
            assert_eq!(ubg.num_extremities(),ubg.num_markers()*2)
        }
        
        for m in ubg.markers() {
            assert!(ubg.degree(tail(m)).is_some());
            assert!(ubg.degree(head(m)).is_some());
        }
        let mset : HashSet<Marker> = ubg.markers().collect();
        for x in ubg.extremities() {
            if !mset.contains(&marker(x)) && x!=TELOMERE {
                panic!("Extremity {x} exists, but its marker does not!");
            }
            if let Some(neighbors) = ubg.adj_neighbors(x) {
                for y in neighbors {
                    if !mset.contains(&marker(y)) && y!=TELOMERE {
                        panic!("Extremity {y} exists, but its marker does not!");
                    }
                let mut has_x = false;
                for z in ubg.adj_neighbors(y).unwrap() {
                    has_x = has_x || (z==x);
                }
                assert!(has_x);
            }
            }
            
        }
    }

    #[test]
    fn test_read_unimog() {
        mtest_read_unimog::<UBG>();
        mtest_read_unimog::<MBG>();
    }
    fn mtest_read_unimog<T>()
    where T : RearrangementGraph
     {
        let ubg =  T::from_unimog("testfiles/test01.ug").expect("File should be readable");
        general_ubg_sanity_check(&ubg);
        for i in 1..=3 {
            assert!(ubg.adj_neighbors(tail(i)).is_some());
            assert!(ubg.adj_neighbors(head(i)).is_some());
        }
        assert!(ubg.adj_neighbors(head(2)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(1),tail(2)])));
        assert!(ubg.adj_neighbors(tail(2)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(3),head(2)])));
        assert!(ubg.adj_neighbors(head(1)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(2)])));
        assert!(ubg.adj_neighbors(head(3)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(1)])));
        assert!(ubg.adj_neighbors(tail(3)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(2)])));
        assert!(ubg.adj_neighbors(tail(1)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(3)])));
    }

    fn collect_neighbors(ubg : &impl RearrangementGraph,xtr : Extremity) -> HashSet<Extremity> {
        ubg.adj_neighbors(xtr).unwrap().collect()
    }

    #[test]
    fn test_unimog_new() {
        mtest_unimog_new::<UBG>();
        mtest_unimog_new::<MBG>();
    }
    fn mtest_unimog_new<T>()
    where T : RearrangementGraph
     {
        let ubg = T::from_unimog("testfiles/test05.ug").expect("File should be readable");
        general_ubg_sanity_check(&ubg);
        for i in ["1","2","5","6","7"] {
            assert!(ubg.name_to_marker(i).is_some());
        }
        for i in ["3","4","+5","+1"] {
            assert!(ubg.name_to_marker(i).is_none());
        }
        let nd1 = ubg.name_to_marker("1").unwrap();
        let nd2 = ubg.name_to_marker("2").unwrap();
        let nd5 = ubg.name_to_marker("5").unwrap();
        let nd6 = ubg.name_to_marker("6").unwrap();
        let nd7 = ubg.name_to_marker("7").unwrap();
        assert!(collect_neighbors(&ubg,tail(nd1)).eq(&HashSet::from([TELOMERE])));
        assert!(collect_neighbors(&ubg,head(nd1)).eq(&HashSet::from([head(nd2),tail(nd6)])));
        assert!(collect_neighbors(&ubg,head(nd2)).eq(&HashSet::from([head(nd1)])));
        assert!(collect_neighbors(&ubg,tail(nd2)).eq(&HashSet::from([tail(nd5)])));
        assert!(collect_neighbors(&ubg,tail(nd5)).eq(&HashSet::from([tail(nd2)])));
        assert!(collect_neighbors(&ubg,head(nd5)).eq(&HashSet::from([TELOMERE])));
        assert!(collect_neighbors(&ubg,tail(nd7)).eq(&HashSet::from([head(nd7)])));
        assert!(collect_neighbors(&ubg,tail(nd6)).eq(&HashSet::from([head(nd1)])));
        assert!(collect_neighbors(&ubg,head(nd6)).eq(&HashSet::from([TELOMERE])));
    }


    #[test]
    fn test_self_cyle() {
        eprintln!("self cycle test ubg");
        mtest_self_cycle::<UBG>();
        eprintln!("self cycle test mbg");
        mtest_self_cycle::<MBG>();
    }

    fn mtest_self_cycle<T>()
    where T : RearrangementGraph
    {
        let mut sizes = HashMap::new();
        let mut adj = HashMap::new();
        let mut nids = HashMap::new();
        nids.insert(String::from("A"), 1);
        sizes.insert(1,1);
        adj.insert(head(1),HashSet::from([head(1)]));
        adj.insert(tail(1), HashSet::from([TELOMERE]));
        adj.insert(TELOMERE, HashSet::from([tail(1)]));
        let g = T::from_hash_maps(sizes, adj, nids);
        let (a,b) = calc_carp_measure_multithread(&g, 1);
        carp_sanity_check(&g, &a.iter().copied().collect(), &b.iter().copied().collect());
        let mut aexp = Vec::new();
        aexp.push((head(1),head(1)));
        assert_eq!(a,aexp);
    }

    #[test]
    fn test_read_gfa_path() {
        let mut mbg = MBG::from_gfa("testfiles/test08.gfa").unwrap();
        general_ubg_sanity_check(&mbg);
        let tels : HashSet<Extremity> = mbg.adj_neighbors(TELOMERE).unwrap().collect();
        let mut expect : HashSet<Extremity> = HashSet::new();
        let m1 = mbg.name_to_marker(&"1").unwrap();
        let m2 = mbg.name_to_marker(&"2").unwrap();
        expect.insert(tail(m1));
        expect.insert(tail(m2));
        assert_eq!(tels,expect);
        mbg.trim_multithread(2, 3);
        general_ubg_sanity_check(&mbg);
        let m1nb : HashSet<Extremity> = mbg.adj_neighbors(head(m1)).unwrap().collect();
        let mut expectnb = HashSet::new();
        expectnb.insert(TELOMERE);
        expectnb.insert(head(m1));
        assert_eq!(m1nb,expectnb);
        mbg.trim_multithread(4, 10);
        assert_eq!(mbg.num_extremities(),0);
        assert_eq!(mbg.adj_neighbors(TELOMERE).map(|x| x.collect()).unwrap_or(HashSet::new()),HashSet::new());
    }


    #[test]
    fn test_read_gfa_walk() {
        let mut mbg = MBG::from_gfa("testfiles/test10.gfa").unwrap();
        general_ubg_sanity_check(&mbg);
        let tels : HashSet<Extremity> = mbg.adj_neighbors(TELOMERE).unwrap().collect();
        let mut expect : HashSet<Extremity> = HashSet::new();
        let m1 = mbg.name_to_marker(&"1").unwrap();
        let m2 = mbg.name_to_marker(&"2").unwrap();
        expect.insert(tail(m1));
        expect.insert(tail(m2));
        assert_eq!(tels,expect);
        mbg.trim_multithread(2, 3);
        general_ubg_sanity_check(&mbg);
        let m1nb : HashSet<Extremity> = mbg.adj_neighbors(head(m1)).unwrap().collect();
        let mut expectnb = HashSet::new();
        expectnb.insert(TELOMERE);
        expectnb.insert(head(m1));
        assert_eq!(m1nb,expectnb);
        mbg.trim_multithread(4, 10);
        assert_eq!(mbg.num_extremities(),0);
        assert_eq!(mbg.adj_neighbors(TELOMERE).map(|x| x.collect()).unwrap_or(HashSet::new()),HashSet::new());
    }


    #[test]
    fn test_read_gfa_path_single() {
        let mbg = MBG::from_gfa("testfiles/test09.gfa").unwrap();
        general_ubg_sanity_check(&mbg);
        let tels : HashSet<Extremity> = mbg.adj_neighbors(TELOMERE).unwrap().collect();
        let mut expect : HashSet<Extremity> = HashSet::new();
        let m1 = mbg.name_to_marker(&"1").unwrap();
        expect.insert(head(m1));
        expect.insert(tail(m1));
        assert_eq!(tels,expect);
    }


    #[test]
    fn test_read_gfa_walk_single() {
        let mbg = MBG::from_gfa("testfiles/test11.gfa").unwrap();
        general_ubg_sanity_check(&mbg);
        let tels : HashSet<Extremity> = mbg.adj_neighbors(TELOMERE).unwrap().collect();
        let mut expect : HashSet<Extremity> = HashSet::new();
        let m1 = mbg.name_to_marker(&"1").unwrap();
        expect.insert(head(m1));
        expect.insert(tail(m1));
        assert_eq!(tels,expect);
    }

    fn scan_sanity_check<T>(g:&T, complexities: &HashMap<Marker,usize>)
    where 
        T : RearrangementGraph
    {
        for m in g.markers() {
            assert!(complexities.contains_key(&m));
        }
    }

    #[test]
    fn test_random_gfa() {
        for i in read_dir("testfiles/random/").expect("W") {
            let x = i.expect("AAAAAAA");
            let tmpval = x.path();
            let gfafile = tmpval.to_str().unwrap();
            eprintln!("Tesing {gfafile}");
            let mut ubg = UBG::from_gfa(gfafile).expect("X");
            ubg.fill_telomeres();
            let mut mbg = MBG::from_gfa(gfafile).expect("Y");
            mbg.fill_telomeres();
            eprintln!("Sanity checking UBG.");
            general_ubg_sanity_check(&ubg);
            eprintln!("Sanity checking MBG.");
            general_ubg_sanity_check(&mbg);
            eprintln!("Eq check");
            equivalence_check(&ubg, &mbg);
            for flt in 0..10 {
                eprintln!("Filtersize or smth: {flt}");
                let mut ubg_ = ubg.clone();
                let mut mbg_ = mbg.clone();
                ubg_.trim_singlethread(flt);
                mbg_.trim_singlethread(flt);
                
                eprintln!("Sanity checking UBG.");
                general_ubg_sanity_check(&ubg_);
                eprintln!("Sanity checking MBG.");
                general_ubg_sanity_check(&mbg_);
                eprintln!("Eq check");
                equivalence_check(&ubg_, &mbg_);
                for i in 3..4 {
                    let mut mbg_mtt = mbg.clone();
                    mbg_mtt.trim_multithread(flt, i);
                    general_ubg_sanity_check(&mbg_mtt);
                    equivalence_check(&mbg_, &mbg_mtt);
                    let mut mbg_mttt = mbg.clone();
                    mbg_mttt.trim_any(flt, i);
                    general_ubg_sanity_check(&mbg_mttt);
                    equivalence_check(&mbg_, &mbg_mttt);
                }
            }
           
        }
    }

    #[test]
    fn test_carp_random_gfa() {
        for i in read_dir("testfiles/random/").expect("W") {
            let x = i.expect("AAAAAAA");
            let tmpval = x.path();
            let gfafile = tmpval.to_str().unwrap();
            eprintln!("Tesing {gfafile}");
            let mut ubg = UBG::from_gfa(gfafile).expect("X");
            ubg.fill_telomeres();
            let mut mbg = MBG::from_gfa(gfafile).expect("Y");
            mbg.fill_telomeres();
            let (c,uc) = calc_carp_measure_naive(&ubg);
            let (cc,ucc) = calc_carp_measure_naive(&mbg);
            let ac = carp_measure_from_adjacencies(&ubg.iter_adjacencies().collect());
            let acc = carp_measure_from_adjacencies(&mbg.iter_adjacencies().collect());
            eprintln!("diff ja lol ey {:?}",c.symmetric_difference(&cc));
            assert_eq!(c,cc);
            assert_eq!(uc,ucc);
            assert_eq!(ac,acc);
            assert_eq!(c.len(),ac);
            assert_eq!(cc.len(),acc);
            carp_sanity_check(&mbg, &c, &uc);
        }
    }

    #[test]
    fn test_scan_random_gfa() {
        for i in read_dir("testfiles/random/").expect("W") {
            let x = i.expect("AAAAAAA");
            let tmpval = x.path();
            let gfafile = tmpval.to_str().unwrap();
            eprintln!("Scan Tesing {gfafile}");
            let mut mbg = MBG::from_gfa(gfafile).expect("Y");
            mbg.fill_telomeres();
            let compl = scan_graph(&mbg, 30);
            let complm = scan_graph_enum_multithread(&mbg, 30, 7);
            eprintln!("Sanity check sthread");
            scan_sanity_check(&mbg, &compl);
            eprintln!("Sanity check mthread");
            scan_sanity_check(&mbg, &complm);
            assert_eq!(compl,complm);
        }
    }

    #[test]
    fn test_read_gfa_fail() {
        mtest_read_gfa_fail::<UBG>();
        mtest_read_gfa_fail::<MBG>();
    }
    fn mtest_read_gfa_fail<T>()
    where T : RearrangementGraph
     {
        for path in ["test01.ug","testfiles/test04.gfa","testfiles/test06.gfa","testfiles/test07.gfa"] {
            let x = T::from_gfa(path);
            assert!(x.is_err());
        } 
    }

    #[test]
    fn test_emoji() {
        mtest_emoji::<UBG>();
        mtest_emoji::<MBG>();
    }
    fn mtest_emoji<T>()
    where T : RearrangementGraph
    {
        let mut ubg = T::from_gfa("testfiles/test03.gfa").expect("File should be readable");
        //eprintln!("nodes {:?}",ubg.node_ids);
        //eprintln!("sizes {:?}",ubg.node_sizes);
        ubg.fill_telomeres();
        general_ubg_sanity_check(&ubg);
        //assert!(ubg.adjacencies.get(&tail(1)).expect(".").eq(&HashSet::from([tail(2)])));
        let h1n : HashSet<Extremity> = ubg.adj_neighbors(tail(1)).unwrap().collect();
        assert!(h1n.eq(&HashSet::from([tail(2)])));
        let t2n : HashSet<Extremity> = ubg.adj_neighbors(tail(2)).unwrap().collect();
        assert!(t2n.eq(&HashSet::from([tail(1)])));
        let ttn : HashSet<Extremity> = ubg.adj_neighbors(TELOMERE).unwrap().collect();
        assert!(ttn.eq(&HashSet::from([head(1),head(2)])));
    }

    #[test]
    fn test_read_gfa() {
        mtest_read_gfa::<UBG>();
        mtest_read_gfa::<MBG>();
    }
    fn mtest_read_gfa<T>()
    where T : RearrangementGraph
     {
        let mut ubg = T::from_gfa("testfiles/test02.gfa").expect("File should be readable");
        ubg.fill_telomeres();
        general_ubg_sanity_check(&ubg);
        for i in 1..=4 {
            assert!(ubg.adj_neighbors(tail(i)).is_some());
            assert!(ubg.adj_neighbors(head(i)).is_some());
        }
        let correct_sizes = HashMap::from([(1,5),(2,2),(3,0),(4,6)]);
        for i in 1..=4 {
            assert_eq!(ubg.node_size(i).unwrap(),*correct_sizes.get(&i).expect("."));
        }
        assert!(ubg.adj_neighbors(head(1)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(4)])));
        assert!(ubg.adj_neighbors(head(2)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(4)])));
        assert!(ubg.adj_neighbors(tail(4)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(1),head(2)])));
        assert!(ubg.adj_neighbors(tail(3)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(2)])));
        assert!(ubg.adj_neighbors(tail(2)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(3)])));
        assert!(ubg.adj_neighbors(head(3)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(4)])));
        assert!(ubg.adj_neighbors(tail(1)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([head(4)])));
        assert!(ubg.adj_neighbors(head(4)).unwrap().collect::<HashSet<_>>().eq(&HashSet::from([tail(1),head(3)])));
    }

    fn carp_multithread_sanity_check(ubg : &impl RearrangementGraph, contested : &HashSet<Adjacency>, uncontested : &HashSet<Adjacency>) {
        for i in 1..8 {
            let (c, uc) = calc_carp_measure_multithread(ubg, i);
            assert_eq!(find_dups(&c),Vec::new());
            assert_eq!(find_dups(&uc),Vec::new());
            eprintln!("El difference: {:?}",contested.symmetric_difference(&c.iter().copied().collect()));
            assert_eq!(c.len(),contested.len());
            assert_eq!(uc.len(),uncontested.len());
            let cs : HashSet<Adjacency> = c.iter().copied().collect();
            let ucs : HashSet<Adjacency> = uc.iter().copied().collect();
            assert_eq!(cs,contested.clone());
            assert_eq!(ucs,uncontested.clone());
        }
    }

    fn carp_sanity_check(ubg : &impl RearrangementGraph, contested : &HashSet<Adjacency>, uncontested : &HashSet<Adjacency>) {
        let mut all = HashSet::new();
        for (x,y) in ubg.iter_adjacencies() {
            all.insert(canonicize((x,y)));
            
        }
        //eprintln!("Mief: {:?}",contested.difference(&all).sorted());
        assert!(contested.is_subset(&all));
        assert!(uncontested.is_subset(&all));
        assert!(contested.is_disjoint(uncontested));
        let union : HashSet<Adjacency> = contested.union(uncontested).cloned().collect();
        //eprintln!("Sym diff: {:?}",all.symmetric_difference(&union));
        assert_eq!(all,union);
        for (x,y) in contested {
            assert!(ubg.adj_neighbors(*x).unwrap().collect::<Vec<_>>().len()>1 
                || ubg.adj_neighbors(*y).unwrap().collect::<Vec<_>>().len()>1 || *x==*y);
            assert!(*x!=TELOMERE && *y!=TELOMERE);
        }
        for (x,y) in uncontested {
            assert!((ubg.adj_neighbors(*x).unwrap().collect::<Vec<_>>().len()==1 
                && ubg.adj_neighbors(*y).unwrap().collect::<Vec<_>>().len()==1 && *x!=*y)
                || *x==TELOMERE || *y==TELOMERE);
            assert!(ubg.adj_neighbors(*x).unwrap().collect::<HashSet<_>>().contains(y)
                && ubg.adj_neighbors(*y).unwrap().collect::<HashSet<_>>().contains(x))
        }
        carp_multithread_sanity_check(ubg, contested, uncontested);

    }

    #[test]
    fn test_carp_eva() {
        let ubg = UBG::from_gfa("testfiles/test02.gfa").expect("File should be readable");
        let mbg = MBG::from_gfa("testfiles/test02.gfa").expect("File should be readable");
        let (contested,uncontested) = calc_carp_measure_naive(&ubg);
        let (mcontested,muncontested) = calc_carp_measure_naive(&mbg);
        assert_eq!(contested,mcontested);
        assert_eq!(uncontested,muncontested);
        carp_sanity_check(&ubg, &contested, &uncontested);
        carp_sanity_check(&mbg, &mcontested, &muncontested);
        assert!(mcontested.len()==4);
        assert!(muncontested.len()==1);
    }

    #[test]
    fn test_back_map() {
        let m = HashMap::from([(1,2),(3,4),(5,6)]);
        let r = reverse_map(&m);
        for (a,b) in m.iter() {
            assert!(r.get(b).expect(".")==a)
        }
        for (a,b) in r.iter() {
            assert!(m.get(b).expect(".")==a)
        }
    }


    fn equivalence_check(g1 : &impl RearrangementGraph, g2 : &impl RearrangementGraph) {
        let a1 : HashSet<Adjacency> = g1.iter_adjacencies().map(|x| canonicize(x)).collect();
        let a2 : HashSet<Adjacency> = g2.iter_adjacencies().map(|x| canonicize(x)).collect();
        let mut mp = g1.marker_names();
        let mp2 = g2.marker_names();
        mp.extend(mp2);
        let a_u1 : HashSet<Adjacency> = a1.difference(&a2).copied().collect();
        eprintln!("{a_u1:?}");
        let a_u2 : HashSet<Adjacency> = a2.difference(&a1).copied().collect();
        eprintln!("{a_u2:?}");
        let a_u : HashSet<String> = a1.symmetric_difference(&a2).copied().map(|a| pretty_adjacency(&mp, a)).collect();
        eprintln!("{a_u:?}");
        assert!(a_u.len()==0);
        let m1 : HashSet<Marker> = g1.markers().collect();
        let m2 : HashSet<Marker> = g2.markers().collect();
        eprintln!("{m1:?}");
        assert_eq!(m1,m2);
        let ms1 : HashSet<(Marker,usize)> = m1.iter().map(|x| 
            {
            //eprintln!("{x}");
            (*x,g1.node_size(*x).unwrap())
            }).collect();
        let ms2 : HashSet<(Marker,usize)> = m2.iter().map(|x| (*x,g2.node_size(*x).unwrap())).collect();
        assert_eq!(ms1,ms2);
    }

    #[test]
    fn test_trim_ypestis() {
        let mut mbg = MBG::from_gfa("testfiles/test_ypestis.gfa").expect("File should be readable");
        let mut ubg = UBG::from_gfa("testfiles/test_ypestis.gfa").expect("File should be readable");
        mbg.fill_telomeres();
        ubg.fill_telomeres();
        equivalence_check(&mbg, &ubg);
        eprintln!("Pre-Testing MBG");
        general_ubg_sanity_check(&mbg);
        eprintln!("Pre-Testing UBG");
        general_ubg_sanity_check(&ubg);
        for i in 32..35 {
            let mut ubg = ubg. clone();
            ubg.trim_singlethread(i);
            for nt in 1..9 {
                let mut mbg = mbg. clone();
                mbg.trim_multithread(i,nt);
                eprintln!("Testing MBG t:{i}");
                general_ubg_sanity_check(&mbg);
                eprintln!("Testing UBG t:{i}");
                general_ubg_sanity_check(&ubg);
                equivalence_check(&mbg, &ubg);
            }
            
        }
    }

    #[test]
    fn test_multithread_trimming_empty() {
        let mut mbg = MBG::from_gfa("testfiles/test_16s.gfa").expect("File should be readable");
        mbg.fill_telomeres();
        mbg.trim_multithread(ONE_MILLION,2);
        general_ubg_sanity_check(&mbg);
    }





    #[test]
    fn test_scan() {
        let mut mbg = MBG::from_gfa("testfiles/test12.gfa").unwrap();
        mbg.fill_telomeres();
        general_ubg_sanity_check(&mbg);
        let m1 = mbg.name_to_marker("1").unwrap();
        let m2 = mbg.name_to_marker("2").unwrap();
        let m3 = mbg.name_to_marker("3").unwrap();
        let m4 = mbg.name_to_marker("4").unwrap();
        let m5 = mbg.name_to_marker("5").unwrap();
        let m6 = mbg.name_to_marker("6").unwrap();
        let m7 = mbg.name_to_marker("7").unwrap();
        let m8 = mbg.name_to_marker("8").unwrap();
        let m9 = mbg.name_to_marker("9").unwrap();
        let m10 = mbg.name_to_marker("10").unwrap();

        //Test for emptiness
        let adj = adjacency_neighborhood(m2,24,&mbg);
        assert_eq!(adj,HashSet::new());


        //Test emptiness for self loops
        let adj = adjacency_neighborhood(m1,24,&mbg);
        assert_eq!(adj,HashSet::new());

        //Test complete component
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m2,500,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m2),tail(m4)),
            (head(m4),head(m2)),
            (head(m1),tail(m2)),
            (tail(m1),tail(m1)),
            (head(m3),tail(m1)),
            (tail(m3),tail(m3)),
            (tail(m3),tail(m2)),
            
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test cutoff
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m2,100,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m2),tail(m4)),
            (head(m4),head(m2)),
            (head(m1),tail(m2)),
            (tail(m3),tail(m2))
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test no parallel
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m5,1500,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m5),tail(m6)),
            (head(m5),head(m10)),
            (tail(m5),tail(m10)),
            (head(m9),tail(m10)),
            (tail(m9),tail(m8)),
            (head(m6),TELOMERE),
            (head(m8),TELOMERE)
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test directed in
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m9,1500,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m5),head(m10)),
            (tail(m5),tail(m10)),
            (head(m9),tail(m10)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE)
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test shortened
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m8,200,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m7),tail(m8)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE)
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);


        //Test shortened
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m8,400,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m7),tail(m8)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE),
            (tail(m7),tail(m6)),
            (head(m9),tail(m10))
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);

        //Test shortened
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m8,650,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m7),tail(m8)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE),
            (tail(m7),tail(m6)),
            (head(m9),tail(m10)),
            (head(m5),head(m10)),
            (head(m6),TELOMERE)
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);

        //Test shortened
        let adj : HashSet<Adjacency> = adjacency_neighborhood(m8,691,&mbg).iter().map(|x| canonicize(*x)).collect();
        let expect : HashSet<Adjacency> = [
            (head(m7),tail(m8)),
            (tail(m9),tail(m8)),
            (head(m8),TELOMERE),
            (tail(m7),tail(m6)),
            (head(m9),tail(m10)),
            (head(m5),head(m10)),
            (head(m6),TELOMERE),
            (tail(5),tail(10))
        ].iter().map(|x| canonicize(*x)).collect();
        assert_eq!(adj,expect);
    }



#[test]
fn test_hisogram() {
    let mut complexities = HashMap::new();
    let hist_empty = histogram(&complexities);
    assert!(hist_empty.is_empty());
    complexities.insert(1, 2);
    complexities.insert(3, 2);
    complexities.insert(4, 3);
    let mut expect = HashMap::new();
    expect.insert(2,2);
    expect.insert(3, 1);
    assert_eq!(expect,histogram(&complexities));
}


#[test]
fn test_percentile() {
    let mut complexities = HashMap::new();
    complexities.insert(1, 2);
    complexities.insert(2,40);
    complexities.insert(3, 2);
    complexities.insert(4, 3);
    complexities.insert(5, 6);
    complexities.insert(6, 1);
    complexities.insert(7, 4);
    complexities.insert(8, 5);
    complexities.insert(9,5);
    complexities.insert(10, 5);
    /*let mut prev = HashSet::new();
    //lower test loop
    for i in 0..10 {
        let tp = top_percentile(&complexities, 0.0, (i as f64/10.0));
        let tps = tp.iter().copied().collect();
        assert!(prev.is_subset(&tps));
        prev = tps;
    }

    let mut prev = HashSet::new();
    //upper test loop
    for i in 1..11 {
        let tp = top_percentile(&complexities, 1.0 - (i as f64/10.0), 1.0);
        let tps = tp.iter().copied().collect();
        assert!(prev.is_subset(&tps));
        prev = tps;
    }*/

    let tps :HashSet<usize> = top_percentile(&complexities, 0.0, 0.3).iter().copied().collect();
    let exp = [6,1,3].iter().copied().collect();
    assert_eq!(tps,exp);

    let tps :HashSet<usize> = top_percentile(&complexities, 0.0, 1.0).iter().copied().collect();
    let exp = (1..11).collect();
    assert_eq!(tps,exp);


    let tps :HashSet<usize> = top_percentile(&complexities, 0.0, 0.9).iter().copied().collect();
    let mut exp : HashSet<usize> = (1..11).collect();
    exp.remove(&2);
    assert_eq!(tps,exp);

    let tps :HashSet<usize> = top_percentile(&complexities, 0.1, 0.2).iter().copied().collect();
    let exp = [1,3].iter().copied().collect();
    assert_eq!(tps,exp);


    let tps :HashSet<usize> = top_percentile(&complexities, 0.2, 0.5).iter().copied().collect();
    let exp = [4,7].iter().copied().collect();
    assert_eq!(tps,exp);

    let tps :HashSet<usize> = top_percentile(&complexities, 0.8, 1.0).iter().copied().collect();
    let exp = [2,5].iter().copied().collect();
    assert_eq!(tps,exp);


    let tps :HashSet<usize> = top_percentile(&complexities, 0.9, 1.0).iter().copied().collect();
    let exp = [2].iter().copied().collect();
    assert_eq!(tps,exp);

}