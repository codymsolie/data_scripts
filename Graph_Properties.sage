# Graph_Properties.sage
#
# this script computes properties of
# the derangement graph associated with
# the group (its Common properties)
#
# *all* of these properties are computed in
# at most O(n) time (linear in number of 
# eigenvalues). Eigenvalues are passed to 
# this class as (value, multiplicity) pairs
#
# *SORTED DECREASING BY VALUE*
#

from collections import defaultdict
import numpy as np

class Graph_Properties:
    def __init__(self, common):
        self.num_vertices = common.order
        self.eigenvalues = common.eigenvalues
        self.max_eigenvalue = common.max_eigenvalue
        self.min_eigenvalue = common.min_eigenvalue
        self.is_union = self._get_is_union()
        self.is_join = self._get_is_join()
        self.is_cmp = self._get_is_cmp()
        self.is_pm_join = self._get_is_pm_join()
        self.is_cograph = self._get_is_cograph()

    # function: _get_is_union
    #
    # return type: boolean
    #
    # the derangement graph is a union if the
    # multiplicity of its largest eigenvalue
    # is greater than one

    def _get_is_union(self):
        return self.eigenvalues[0][1] > 1

    # function: _get_is_join
    #
    # return type: boolean
    #
    # the derangement graph is a join if 
    # (max_evalue - min_evalue) == num_vertices

    def _get_is_join(self):
        return (self.max_eigenvalue - self.min_eigenvalue) == self.num_vertices

    # function: _get_is_cmp
    #
    # return type: boolean
    #
    # the derangement graph is a complete 
    # multipartite graph if it is a join,
    # and there are exactly three distinct
    # eigenvalues (omitting mult.), namely:
    #
    # {(degree, m1), (0, m2), (\tau, m3)}

    def _get_is_cmp(self):
        eigenvalues = set(map(itemgetter(0), self.eigenvalues)) # removes multiplicity
        desired = set([int(self.max_eigenvalue), int(0), int(self.min_eigenvalue)])
        return (self.is_join and eigenvalues == desired)

    # function: _get_is_pm_join
    #
    # return type: boolean
    #
    # the derangement graph is a join
    # of perfect matchings if it is 
    # a join, and the eigenvalues are
    # as follows:
    #
    #   {-2k+1, -1, 1, 2kl-2k+1}

    def _get_is_pm_join(self):
      if self.is_join:
        evalues = sorted([pair[0] for pair in self.eigenvalues])
        k = (1 - self.min_eigenvalue) / 2
        l = self.eigenvalues[0][1] + 1
        desired = [-(2*k)+1, -1, 1, (2*k*l)-(2*k)+1]
        if desired == evalues:
            return True
      return False

    # function: _get_is_cograph
    #
    # return type: boolean
    #
    # the derangement graph is a cograph
    # if it can be written as consecutive
    # joins and unions of K_1, the trivial
    # graph with a single vertex.
    #
    # if the graph is a union, we can 'reduce'
    # it by dividing off each multiplicity 
    # until the multiplicity of the degree is 
    # equal to one
    #
    # if the graph is a join, we can 'reduce'
    # it by taking its complement, which must
    # be a disjoint union, reducing that as
    # above, then flipping back out of the
    # complement. The eigenvalues are kept
    # intact by a theorem on regular graphs, 
    # and the eigenvalues of their complements.

    def _get_is_cograph(self):

        # create new instance to preserve class data
        eigenvalues = self.eigenvalues
        # this next line changes (immutable) pairs
        # to (mutable) lists
        eigenvalues = list(map(list, eigenvalues))
        max = self.max_eigenvalue
        min = self.min_eigenvalue
        v = self.num_vertices

        def get_complement_evalues(eigenvalues, v):
            # this is the v-d-1 case
            eigenvalues[0][0] = v-eigenvalues[0][0]-1

            # these are the -1-l cases for complement
            for pair in eigenvalues[1:]:
                pair[0] = -1-pair[0]

            # consolidate multiplicities (the right way????)
            dict = defaultdict(int)
            for e, m in eigenvalues:
                dict[e]+=m

            # sort the eigenvalues decreasing by first
            # in pair
            new_evalues = list(dict.items())
            # convert tuple pairs (immutable) to lists
            new_evalues = np.array(new_evalues).tolist()

            new_evalues.sort(reverse = True, key=lambda x: x[0])

            return new_evalues

        def computation(eigenvalues, max, min, v):
            # union case
            if eigenvalues[0][1] > 1:
                num_components = eigenvalues[0][1]
                union_part_evalues = list(map(lambda pr: [pr[0], int(pr[1]/num_components)], eigenvalues))
                return computation(union_part_evalues, max, min, v/num_components)
            # join case
            elif (max - min) == v:
                new_degree = v-eigenvalues[0][0]-1
                comp_evalues = get_complement_evalues(eigenvalues, v)
                num_components = comp_evalues[0][1]
                # reduce multiplicities in complement (union)
                lower_comp_evalues = list(map(lambda pr: [pr[0], int(pr[1]/num_components)], comp_evalues))
                join_part_evalues = get_complement_evalues(lower_comp_evalues, v/num_components)
                new_max = join_part_evalues[0][0]
                new_min = join_part_evalues[-1][0]
                return computation(join_part_evalues, new_max, new_min, v/num_components)
            # irreducible case
            else:
                return v == 1

        return computation(eigenvalues, max, min, v)
