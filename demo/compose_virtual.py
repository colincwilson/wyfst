# Demo composition of wfst1 and virtual wfst2.
import os, sys

from wyfst import *
from pynini import Weight

# Machine wfst1
wfst1 = wywrapfst.accep('a b c', isymbols=None)
print(wfst1.info())


# Virtual machine wfst2
# States represent optional FIFO buffered inputs.
def wfst2_func(src2, t1_olabel, t1_weight):
    (q, buff) = src2
    ret = []
    if t1_olabel == config.epsilon:
        if q == 'I':
            ret.append((('I', buff), t1_olabel, t1_weight))
            ret.append((('B', buff), t1_olabel, t1_weight))
        if q == 'B':
            ret.append((('B', buff), t1_olabel, t1_weight))
            ret.append((('H', buff), t1_olabel, t1_weight))
    elif q == 'B':
        buff_ = buff + ' ' + t1_olabel
        ret.append(((q, buff_), t1_olabel, t1_weight))
    elif q in ['I', 'H']:
        ret.append((src2, t1_olabel, t1_weight))
    return ret


initial2 = ('I', '')


def final2_func(src2):
    return Weight.one('tropical')


# Composition output.
wfst = wywrapfst.compose_virtual(wfst1,
                                 wfst2_func,
                                 initial2,
                                 final2_func,
                                 verbose=True)
wfst.info()
wfst.draw('fig/compose_virtual.dot')
