import unittest
import utils as u
import scj_carp as c
import random as r

def xtr(y):
    if y=='t':
        return u.EXTREMITY_TAIL
    elif y=='h':
        return u.EXTREMITY_HEAD
    else:
        return u.TELO

def hr2e(x):
    xs = x.split(',')
    a=(xs[0][0:-1],xtr(xs[0][-1]))
    b=(xs[1][0:-1],xtr(xs[1][-1]))
    return tuple(sorted([a,b]))


TLNODE = (u.TELO,u.TELO)

#linear graph
def graph1():
        l = []
        g1 = ('>G',[(u.CHR_LINEAR,[(u.ORIENT_POSITIVE,'1'),(u.ORIENT_POSITIVE,'2'),(u.ORIENT_NEGATIVE,'3')])])
        g2 = ('>G',[(u.CHR_LINEAR,[(u.ORIENT_POSITIVE,'3'),(u.ORIENT_POSITIVE,'2'),(u.ORIENT_POSITIVE,'2')])])
        l.append(g1)
        l.append(g2)
        mbpg = u.build_multi_bp_graph(l)
        return mbpg

#circular graph
def graph2():
        g1 = ('>G',[(u.CHR_CIRCULAR,[(u.ORIENT_POSITIVE,'1'),(u.ORIENT_POSITIVE,'2'),(u.ORIENT_NEGATIVE,'2')])])
        g2 = ('>G',[(u.CHR_CIRCULAR,[(u.ORIENT_POSITIVE,'3'),(u.ORIENT_POSITIVE,'2')])])
        g3 = ('>G',[(u.CHR_CIRCULAR,[(u.ORIENT_NEGATIVE,'6')])])
        l=[g1,g2,g3]
        mbpg = u.build_multi_bp_graph(l)
        return mbpg

#mixed graph
def graph3():
        g1 = ('>G',[(u.CHR_CIRCULAR,[(u.ORIENT_POSITIVE,'1'),(u.ORIENT_POSITIVE,'2'),(u.ORIENT_NEGATIVE,'2')])])
        g2 = ('>G',[(u.CHR_CIRCULAR,[(u.ORIENT_POSITIVE,'3'),(u.ORIENT_POSITIVE,'2')])])
        g3 = ('>G',[(u.CHR_LINEAR,[(u.ORIENT_NEGATIVE,'1')])])
        l=[g1,g2,g3]
        mbpg = u.build_multi_bp_graph(l)
        return mbpg



class TestMBPG(unittest.TestCase):
    def test_graph_linear(self):
        mbpg = graph1()
        nodes = set(mbpg.nodes())
        edges = set([tuple(sorted(e)) for e in mbpg.edges()])
        self.assertEqual(nodes,set([(x,y) for x in "123" for y in [u.EXTREMITY_HEAD,u.EXTREMITY_TAIL]]+[(u.TELO,u.TELO)]))
        self.assertEqual(edges,set([hr2e(x) for x in ["oo,1t","1h,2t","2h,3h","3t,oo","oo,3t","3h,2t","2h,2t","2h,oo"]]))
    def test_graph_circular(self):
        mbpg = graph2()
        nodes = set(mbpg.nodes())
        edges = set([tuple(sorted(e)) for e in mbpg.edges()])
        self.assertEqual(nodes,set([(x,y) for x in "1236" for y in [u.EXTREMITY_HEAD,u.EXTREMITY_TAIL]]))
        self.assertEqual(edges,set([hr2e(x) for x in ["1h,2t","2h,2h","2t,1t","3h,2t","2h,3t","6h,6t"]]))
    def test_graph_mixed(self):
        mbpg = graph3()
        nodes = set(mbpg.nodes())
        edges = set([tuple(sorted(e)) for e in mbpg.edges()])
        self.assertEqual(nodes,set([(x,y) for x in "123" for y in [u.EXTREMITY_HEAD,u.EXTREMITY_TAIL]]+[(u.TELO,u.TELO)]))
        self.assertEqual(edges,set([hr2e(x) for x in ["1h,2t","2h,2h","2t,1t","3h,2t","2h,3t","1t,oo","1h,oo"]]))


class TestCARP(unittest.TestCase):
    def test_graph_linear(self):
        lmbpg = graph1()
        m, con,unc = c.calc_carp_index(lmbpg,get_edge_partition=True)
        self.sanity_checks(lmbpg, m, con, unc)
        self.assertEqual(m,4)
        cons = set([tuple(sorted(e)) for e in con])
        self.assertEqual(cons,set([hr2e(x) for x in ["1h,2t","2h,3h","3h,2t","2h,2t"]]))

    def test_graph_circular(self):
        cmbpg = graph2()
        #print(cmbpg.nodes())
        #print(cmbpg.edges())
        m, con,unc = c.calc_carp_index(cmbpg,get_edge_partition=True)
        self.sanity_checks(cmbpg,m,con,unc)
        cons = set([tuple(sorted(e)) for e in con])
        self.assertEqual(cons,set([hr2e(x) for x in ["1h,2t","2h,2h","2t,1t","3h,2t","2h,3t"]]))
    def test_graph_mixed(self):
        lmbpg = graph3()
        m, con,unc = c.calc_carp_index(lmbpg,get_edge_partition=True)
        self.sanity_checks(lmbpg, m, con, unc)
        self.assertEqual(m,5)
        cons = set([tuple(sorted(e)) for e in con])
        self.assertEqual(cons,set([hr2e(x) for x in ["1h,2t","2h,2h","2t,1t","3h,2t","2h,3t"]]))
    def test_identity(self):
         chrsa = [(u.CHR_LINEAR,[(u.ORIENT_NEGATIVE,'1'),(u.ORIENT_POSITIVE,'4'),(u.ORIENT_POSITIVE,'6'),]),(u.CHR_CIRCULAR,[(u.ORIENT_POSITIVE,'3')])]
         chrsb = [(u.CHR_LINEAR,[(u.ORIENT_POSITIVE,'7')])]
         gnms = []
         for i,nm in enumerate("ABCDE"):
             if i%2==0:
                gnms.append((nm,chrsa))
             else:
                 gnms.append((nm,chrsb))
         mbpg = u.build_multi_bp_graph(gnms)
         m,con,unc = c.calc_carp_index(mbpg,get_edge_partition=True)
         self.sanity_checks(mbpg,m,con,unc)
         self.assertEqual(m,0)
    def test_oneloop(self):
        g1 = ('>G',[(u.CHR_LINEAR,[(u.ORIENT_POSITIVE,'1'),(u.ORIENT_NEGATIVE,'1')])])
        mbpg = u.build_multi_bp_graph([g1])
        m,con,unc = c.calc_carp_index(mbpg,get_edge_partition=True)
        self.sanity_checks(mbpg,m,con,unc)
        self.assertEqual(m,1)
        self.assertEqual(con,[hr2e("1h,1h")])
    
    def test_random(self):
        r.seed(42)
        REPLICATES = 1000
        for _ in range(REPLICATES):
            ngenomes = r.randint(1,5)
            nfams = r.randint(10,100000)
            gnms = []
            for _ in range(ngenomes):
                nchroms = r.randint(1,20)
                gnmname=''.join(r.choices("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"))
                chmrms = []
                for _ in range(nchroms):
                    nmarkers = r.randint(1,1000)
                    chrtype = r.choice([u.CHR_CIRCULAR,u.CHR_LINEAR])
                    chrstr =[]
                    for _ in range(nmarkers):
                        chrstr.append((r.choice([u.ORIENT_NEGATIVE,u.ORIENT_POSITIVE]),str(r.randint(1,nfams))))
                    chmrms.append((chrtype,chrstr))
                gnms.append((gnmname,chmrms))
        mbpg = u.build_multi_bp_graph(gnms)
        m,con,unc = c.calc_carp_index(mbpg,get_edge_partition=True)
        cons = set([tuple(sorted(e)) for e in con])
        uncs = set([tuple(sorted(e)) for e in unc])
        self.sanity_checks(mbpg,m,con,unc)
        gg = r.choice(gnms)
        mbpg_ = u.build_multi_bp_graph(gnms+[gg])
        m_,con_,unc_ = c.calc_carp_index(mbpg_,get_edge_partition=True)
        cons_ = set([tuple(sorted(e)) for e in con_])
        uncs_ = set([tuple(sorted(e)) for e in unc_])
        self.assertEqual(m,m_)
        self.assertEqual(cons,cons_)
        self.assertEqual(uncs,uncs_)


    def sanity_checks(self, mbpg, m, con, unc):
        cons = set([tuple(sorted(e)) for e in con])
        uncs = set([tuple(sorted(e)) for e in unc])
        #print(cons)
        #print(uncs)
        #sanity checks
        self.assertEqual(len(cons.intersection(uncs)),0)
        alle = cons.union(uncs)
        edges = set([tuple(sorted(e)) for e in mbpg.edges()])
        self.assertEqual(alle,edges)
        self.assertEqual(m,len(cons))
        for e in uncs:
            self.assertNotEqual(e[0],e[1])
            v,w = e
            if v != TLNODE and w!= TLNODE:
                self.assertEqual(mbpg.degree[v],1)
                self.assertEqual(mbpg.degree[w],1)
        for e in cons:
             self.assertNotIn(u.TELO,e)
             v,w = e

             self.assertNotEqual((mbpg.degree[v],mbpg.degree[w]),(1,1))
        
        



if __name__ == '__main__':
    unittest.main()