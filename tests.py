import unittest
import utils as u
import scj_carp as c


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


TLNODE = hr2e('oo')

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
        mbpg = graph1()
        m, con,unc = c.calc_carp_index(mbpg,get_edge_partition=True)
        cons = set([tuple(sorted(e)) for e in con])
        uncs = set([tuple(sorted(e)) for e in unc])
        #sanity checks
        self.assertEqual(len(cons.intersection(uncs)),0)
        alle = cons.union(uncs)
        edges = set([tuple(sorted(e)) for e in mbpg.edges()])
        for e in edges:
            if TLNODE in e:
                continue
            self.assertIn(e,alle)
        for e in alle:
             self.assertIn(e,edges)
        



if __name__ == '__main__':
    unittest.main()