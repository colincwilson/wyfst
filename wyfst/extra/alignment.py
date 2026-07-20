# String alignment weighted by feature similarity.
import re, sys

import panphon
from panphon.distance import Distance
#from panphon.distance import unweighted_substitution_cost

from phonopy.str_util import str_split
import wyfst
from wyfst import config as wyconfig
from wyfst import Weight


class Aligner:

    def __init__(self, sigma, insertion_cost=0.5, deletion_cost=0.5):
        symbols, _ = wyconfig.init({'sigma': sigma})
        self.symbols = symbols
        self.insertion_cost = insertion_cost
        self.deletion_cost = deletion_cost
        self.seg_distance = seg_distance = Distance()

        self.A = A = wyfst.empty_transducer(isymbols=symbols)

        q = 1
        for x in sigma:
            ftrs_x = self.seg_distance.fm.word_to_vector_list( \
                x, numeric=True)[0]
            for y in sigma:
                ftrs_y = self.seg_distance.fm.word_to_vector_list( \
                    y, numeric=True)[0]
                cost_xy = seg_distance.unweighted_substitution_cost(
                    ftrs_x, ftrs_y)
                A.add_arc(src=q, ilabel=x, olabel=y, weight=cost_xy, dest=q)

        epsilon = wyconfig.epsilon
        for x in sigma:
            A.add_arc(src=q,
                      ilabel=x,
                      olabel=epsilon,
                      weight=self.deletion_cost,
                      dest=q)
            A.add_arc(src=q,
                      ilabel=epsilon,
                      olabel=x,
                      weight=self.insertion_cost,
                      dest=q)

    def align(self, s1, s2):
        I = wyfst.acceptor(s1, isymbols=self.symbols, add_delim=True)
        O = wyfst.acceptor(s2, isymbols=self.symbols, add_delim=True)
        M1 = wyfst.compose(I, self.A)
        M2 = wyfst.compose(self.A, O)
        M = wyfst.compose(M1, M2)
        align_best = wyfst.shortestpath(M)
        ret = list(align_best.iostrings())
        return ret


if __name__ == "__main__":
    s1 = str_split("tristi")
    s2 = str_split("trumisti")

    sigma = s2.split(" ")
    aligner = Aligner(sigma=sigma)

    a = aligner.align(s1, s2)
    print(a)
