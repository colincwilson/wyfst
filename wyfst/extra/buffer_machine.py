# Demo composition of wfst1 and virtual wfst2.
import os, sys

from wyfst import config as wyconfig
from wyfst.wywrapfst import *
from pynini import Weight

#print(vars(config))
#sys.exit(0)

epsilon = wyconfig.epsilon

weight_one = Weight.one('tropical')  # 0
weight = Weight('tropical', 1.0)  # 1
weight_zero = Weight.zero('tropical')  # +Inf

# Machine wfst1
# inpt = 'w o f e w o'  # Test stim inMarcus et al. 1999
inpt = 'a ng g a ng a ng g a ng'  # Listed in MALINDO as R-penuh
wfst1 = accep(inpt, isymbols=None)
print(wfst1.info())


# Virtual machine wfst2
# Second component of state is optional FIFO buffered input.
def wfst2_func(src2, t1_olabel, t1_weight):
    (q, buff) = src2
    arcs = []
    # Change state on epsilon transitions.
    if t1_olabel == epsilon:
        if q == 'I':
            arcs.append((src2, '<', weight_one, ('B0', buff)))
            # note: distinction between B0 and B filters out empty-buffer paths
        if q == 'B':
            arcs.append((src2, '>', weight_one, ('H', buff)))
        if q == 'H':
            arcs.append((src2, '<', weight_one, ('M', buff)))
        if q == 'M':
            arcs.append((src2, '>', weight_one, ('F', buff)))
    # Scan or buffer hold.
    elif q in ['I', 'H', 'F']:
        arcs.append((src2, t1_olabel, weight, src2))
    # Buffer store.
    elif q in ['B0', 'B']:
        buff_ = t1_olabel if buff == '' else (buff + ' ' + t1_olabel)
        arcs.append((src2, t1_olabel, weight_one, ('B', buff_)))
    # Buffer match.
    elif q == 'M' and buff != '':
        try:
            head, _buff = buff.split(' ', maxsplit=1)
        except:
            head, _buff = buff, ''
        if head == t1_olabel:
            arcs.append((src2, t1_olabel, weight_one, (q, _buff)))
    return arcs


initial2 = ('I', '')


def final2_func(src2):
    (q, buff) = src2
    if q == 'F' and buff == '':
        return Weight.one('tropical')
    return Weight.zero('tropical')


# Composition output.
wfst = compose_virtual( \
    wfst1, wfst2_func, initial2, final2_func, verbose=False)

print('=' * 5 + 'all paths' + '=' * 5)
outs = wfst.ostrings()
for x in outs:
    print(x)

wfst.info()

wfst.relabel_states()
wfst.draw('../../demo/fig/buffer_machine.dot')

# Best path (longest self-match).
print('=' * 5 + 'best path' + '=' * 5)
wfst_best = wfst.shortestpath()
outs_best = wfst_best.ostrings()
for x in outs_best:
    print(x)
