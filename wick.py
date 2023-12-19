#!/bin/python3
import sys
import copy
import itertools
from IPython.display import display, Latex

def do_wick(string):
    debug = False
    
    if debug: print('Input string:',string)

    idx_orb = 0
    idx_op  = 2

    pos = [i for i in range(len(string))]
    # Extraction of the operators
    ## By class of orb
    ### i
    i_op, pos_i_op = extract_op(string,'i',idx_orb,pos)
    if debug: print('i op:',i_op,'\nPos i op:',pos_i_op)
    ### a
    a_op, pos_a_op = extract_op(string,'a',idx_orb,pos)
    if debug: print('a op:',a_op,'\nPos a op:',pos_a_op)
    
    ## By class of op
    ### i+
    i_crea, pos_i_crea = extract_op(i_op,'+',idx_op,pos_i_op)
    if debug: print('i+:',i_crea,'\nPos i+:',pos_i_crea)
    ### i-
    i_anni, pos_i_anni = extract_op(i_op,'-',idx_op,pos_i_op)
    if debug: print('i-:',i_anni,'\nPos i-:',pos_i_anni)
    ### a+
    a_crea, pos_a_crea = extract_op(a_op,'+',idx_op,pos_a_op)
    if debug: print('a+:',a_crea,'\nPos a+:',pos_a_crea)
    ### a-
    a_anni, pos_a_anni = extract_op(a_op,'-',idx_op,pos_a_op)
    if debug: print('a-:',a_anni,'\nPos a-:',pos_a_anni)

    # Build delta
    ## i+ i-
    delta_i = build_delta(i_crea,i_anni,pos_i_crea,pos_i_anni)
    if debug: print('Delta i:', delta_i)
    ## a- a+
    delta_a = build_delta(a_anni,a_crea,pos_a_anni,pos_a_crea)
    if debug: print('Delta a:', delta_a)
    #sys.exit()

    # Creation of all the products of hole deltas
    res_i = combine_one_class_delta(delta_i,i_op,string)
    res_i.insert(0,[[]])
    if debug: print('res_i (possible products of i deltas):',res_i)
    # Creation of all the products of particle deltas
    res_a = combine_one_class_delta(delta_a,a_op,string)
    res_a.insert(0,[[]])
    if debug: print('res_a (possible products of a deltas):',res_a)

    # Combination of the products of kronecker deltas coming from different classes of orbitals
    res = []
    # Iteration over the order (order <=> number of deltas) for i
    for order_i in res_i:
        # Iteration over the list of product of deltas i of a given order
        for list_delta_i in order_i:
            # As for i
            for order_a in res_a:
                # As for i
                for list_delta_a in order_a:
                    # If there is at least one term for a
                    if len(list_delta_a) > 0:
                        # the different deltas a are append after the i ones
                        tmp = copy.deepcopy(list_delta_i)
                        for delta_a in list_delta_a:
                            tmp.append(delta_a)
                        res.append(tmp)
                    # If there is no a delta, we just keep the i ones
                    else:
                        tmp = copy.deepcopy(list_delta_i)
                        res.append(tmp)

    # Ordering of the terms by increasing number of deltas
    max_len = 0
    for elem in res:
        if len(elem) > max_len:
            max_len = len(elem)

    # Reordering of the list of deltas depending on the number of deltas
    tmp_ordered = [[] for i in range(max_len+1)]
    for elem in res:
        tmp_ordered[len(elem)].append(elem)
    
    # Finally the right order
    list_deltas = []
    for l_deltas in tmp_ordered:
        for elem in l_deltas:
            list_deltas.append(elem)
    
    list_signs = []
    for elem in list_deltas:
        list_signs.append(extract_sign(elem,string))

    # copy of the initial string
    list_str = [string for i in range(len(list_signs))]

    # Removing the op of the string implied in a contraction
    acc = []
    for s,d in zip(list_str,list_deltas):
        tmp = []
        for elem in s:
            is_contracted = False
            for delta in d:
                for op in delta:
                    if op == elem:
                        is_contracted = True
                        break
            if (not is_contracted):
                tmp.append(elem)

        acc.append(tmp)
        
    list_str = acc

    if debug:
        for sign,d,s in zip(list_signs,list_deltas,list_str):
            print(sign,d,s)

    return list_signs, list_deltas, list_str

# The idea is to build the string by starting from the list of delta
# From this list we can build a list of pairs of deltas that do not share the same idx
# From this list of pairs of delta we can built a list of triplet of delta
# and so on ...
def combine_one_class_delta(delta_i,i_op,string):
    
    order_delta_i = []
    if len(delta_i) == 0:
        return order_delta_i
        
    acc = []
    if len(delta_i) != 0:
        for d in delta_i:
            acc.append([d])
    else:
        acc.append([(0,0)])
    order_delta_i.append(acc)
   
    for i in range(1,len(i_op)//2+1):
        # List of term
        acc1 = []
        # a term = 1 or many delta that are multiplied
        list_term = order_delta_i[i-1]
        # list of delta of each term, 1 or many delta that are multiplied
        acc2 = []
        for list_delta in list_term:
            list_op = []
            # a single delta
            for delta in list_delta:
               # the op that are contracted with the delta
               for op in delta:
                   list_op.append(op)
                   
            # add delta, one on the existing delta and we look if we can do the contraction in
            # addition do the contractions already done
            for add_delta in delta_i:
                acc3 = copy.deepcopy(list_delta)
                idx_last = find_idx_elem(acc3[-1] ,delta_i)
                idx_add  = find_idx_elem(add_delta,delta_i)
                if idx_add <= idx_last:
                    continue
                is_in = False
                # check if the operator in the delta we want to add is already in another contraction
                for op1 in add_delta:
                    for op2 in list_op:
                        if op1 == op2:
                            is_in = True
                    if is_in: continue
                if is_in: continue
                
                acc3.append(add_delta)
                acc2.append(acc3)
                
        #if there is no i-multiple delta        
        if len(acc2) == 0:
            break
        order_delta_i.append(acc2)

    # sign
    #for list_term in order_delta_i:
    #    for list_delta in list_term:
    #        sign = extract_sign(list_delta,string)

    return order_delta_i

# Product of two lists
# To compute the sign based on a list of delta and the original string of operators            
def extract_sign(list_delta,string):
    # position of the op in the different deltas
    list_pos = []
    for delta in list_delta:
        tmp = []
        for op in delta:
            pos = find_idx_elem(op,string)
            tmp.append(pos)
        list_pos.append((tmp[0],tmp[1]))

    # Number of crossing lines in the contractions
    nb_cross = 0
    for j in range(0,len(list_pos)-1):
        for k in range(j+1,len(list_pos)):
            pi_j = list_pos[j][0]
            pf_j = list_pos[j][1]
            pi_k = list_pos[k][0]
            pf_k = list_pos[k][1]
            # Crossing : ...pi_j ... pi_k ... pf_j ... pf_k...
            if (pf_k > pf_j) and (pi_k < pf_j) and (pi_k > pi_j) :
                nb_cross = nb_cross + 1

    # Number of permutation required to do the contraction
    nb_perm = 0
    for pos in list_pos:
        pi = pos[0]
        pf = pos[1]
        nb_perm = nb_perm + pf-pi-1

    # Final sign
    sign = (-1)**nb_cross * (-1)**nb_perm
    
    return sign

# extract a type of operator based on a pattern at the idxth position in string[:]
def extract_op(string,pattern,idx,pos):
    debug = False
    
    # Check type
    if type(string) != type(['a','b']):
        print('Type mismatch function extract_op arg 1')
        sys.exit()
    if type(pattern) != type('a'):
        print('Type mismatch function extract_op arg 2')
        sys.exit()
    if type(idx) != type(1):
        print('Type mismatch function extract_op arg 3')
        sys.exit()
    if type(pos) != type([1,2]):
        print('Type mismatch function extract_op arg 4')
        sys.exit()

    # Debug
    if debug: print('string:',string)
    if debug: print('pattern:',pattern)
    if debug: print('idx:',idx)
    
    res = []
    new_pos = []
    i = 0
    for elem in string:
        if elem[idx] == pattern:
            res.append(elem)
            new_pos.append(pos[i])
        i = i + 1

    return res, new_pos

# Build all the possible kronecker delta using 2 list of operators
# and their position in the original string of operator
def build_delta(list_op1,list_op2,list_pos1,list_pos2):
    debug = False

    idx_spin = 3
    idx_act  = 4

    if debug: print('List op 1:',list_op1,'\List op 2:',list_op2)
    if debug: print('List pos 1:',list_pos1,'\List pos 2:',list_pos2)

    nb_idx = len(list_op1[0])
    if nb_idx < 4 or nb_idx > 5:
        print('The operators must have at least 4 indexes and maximum 5 indexes.')
        sys.exit()
    for elem in list_op1:
        if len(elem) != nb_idx:
            print('All the operators must share the same number of indexes.')
            sys.exit()

    a1 = ''
    a2 = ' '
    res = []
    for op1, pos1 in zip(list_op1,list_pos1):
        s1 = op1[idx_spin]
        if nb_idx == 5:
            a1 = op1[idx_act]
        for op2, pos2 in zip(list_op2,list_pos2):
            s2 = op2[idx_spin]
            if nb_idx == 5:
                a2 = op2[idx_act]
            # if alpha-beta spin
            if (s1 == 'a' and s2 == 'b') or (s1 == 'b' and s2 == 'a'):
                continue
            # if active-active contraction or not active-not active contraction
            if a1 == a2 and nb_idx > idx_act:
                continue
            if pos2 > pos1:
                res.append([op1,op2])

    if debug: print('Res:',res)
    
    return res

# To search the index of an element in a list
def find_idx_elem(elem,list_elem):
    i = 0
    for d in list_elem:
        if d == elem:
            break
        else:
            i = i + 1
            
    # check
    if i == len(list_elem):
        print('elem not found in find_idx_elem')
        sys.exit()

    return i

# To put creation operator on the right
def put_crea_to_left(sign,string):
    acc = []
    acc_pos = []
    tmp = []
    idx = 2
    string_pos = [i for i in range(len(string))]
    for elem,pos in zip(string,string_pos):
        if elem[idx] == '+':
            acc.append(elem)
            acc_pos.append(pos)
        else:
            tmp.append(elem)

    order = [i for i in range(len(acc))]
    d = 0
    for pi,pf in zip(acc_pos,order):
        d = d + abs(pi-pf)

    sign = pow(-1,d) * sign
    
    for elem in tmp:
        acc.append(elem)
        
    return sign, acc

class Wicked_str():
    def __init__(self,sign,deltas,ops):
        self.sign = sign
        self.deltas = deltas
        self.ops = ops
        self.tex = self.to_latex()

    def to_latex(self):
        sign = self.sign
        deltas = self.deltas
        ops = self.ops
        if sign > 0:
            tex = '+ '
        else:
            tex = '- '

        tx = deltas_to_tex(deltas)
        tex = tex + tx

        if len(ops) > 0:
            tex = tex + '\\left\\{'
            for op in ops:
                #print(op)
                o = latexify(op)
                tex = tex + o
            tex = tex + '\\right\\}_N'
        return tex
        
    def crea_to_left(self):
        self.sign, self.ops = put_crea_to_left(self.sign,self.ops)
        self.tex = self.to_latex()

    def tex_show(self):
        print(self.tex)
        
    def eq_show(self):
        display(Latex(f'${self.tex}$'))

def deltas_to_tex(deltas):
    tex = ''
    for delta in deltas:
        d1 = str(delta[0][1])+ '_{$'+ delta[0][3] + '}'
        d2 = str(delta[1][1])+ '_{$'+ delta[1][3] + '}'
        tex = tex + '\delta('+d1+','+d2+') \ '
        tex = tex.replace('$a','\\alpha')
        tex = tex.replace('$b','\\beta')
        tex = tex.replace('$g','')
    return tex
    
def latexify(op):
    tex = op[0]+'^{'+op[2]+'}'+'_{'+op[1]+'_{$'+op[3]+'}'+'}'
    tex = tex.replace('+','\dagger')
    tex = tex.replace('-','')
    tex = tex.replace('$a','\\alpha')
    tex = tex.replace('$b','\\beta')
    tex = tex.replace('$g','g')

    return tex

# 1: orbital class, i for occupied, a for unoccupied (for Fermi vacuum).
# For the true vacuum, use a.
# 2: orbital label
# 3: operator type, + for creation, - for annihilation
# 4: spin, a for $\alpha$, b for $\beta$, g for general (could be $\alpha$ or $\beta$)
# 5: optional, to avoid contraction between some operators
# (two operators with the same 5th index cannot be contracted together).

if __name__ == "__main__":
    s = ['ip+g','aq-g','ir-g','is+g','at+g','iu-g']
    #s = ['ix+g','iy+g','aw-g','av-g','iq+g','ip-g']

    #From this, we can call the function "do_wick" on s to apply
    # Wick's theorem and generate 3 lists, one for the signs, one
    # for the kronecker delta and one for the normal ordered string
    # containing the remaining uncontracted operators (WARNING:
    # the operators are not put in normal order in these strings,
    # but you can reorder them since it will only change the sign.
    # That's why we write them as $\{...\}_N$).
    list_sign, list_deltas, list_string = do_wick(s)

    for sign,deltas,string in zip(list_sign, list_deltas, list_string):
        # Creates an object Wicked_str to print or display latex code 
        obj = Wicked_str(sign,deltas,string)
        # to print the each element of the result with latex format
        obj.tex_show()
        # to show the latex equation in a Jupyter Notebook
        #obj.eq_show()
