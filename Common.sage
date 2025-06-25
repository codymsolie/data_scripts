# Common.sage
# 
# collects general properties of the 
# provided group
#
# called upon by Generator, with
# instances passed to EKR_Determiner
# and Intersection_Density (.sage)

from itertools import groupby
from operator import itemgetter

class Common: 
    def __init__(self, group):
        self.group = group
        self.degree = group.degree()
        self.order = group.order()
        self.identity = group.one()
        self.transitivity = gap.Transitivity(group)
        self.conjugacy_classes = group.conjugacy_classes()
        self.subgroups = group.conjugacy_classes_subgroups()
        self.characters = group.irreducible_characters()
        self.size_of_stabilizer = self.group.stabilizer(1).order()
        self.derangement_classes = self._get_derangement_classes()
        self.eigenvalues = self._get_eigenvalues()
        self.max_eigenvalue = self._get_max_eigenvalue()
        self.min_eigenvalue = self._get_min_eigenvalue()
        self.n_cliques = self._get_n_cliques()
        self.stabilizer_sized_cocliques = self._get_stabilizer_sized_cocliques()
        self.larger_than_stabilizer_cocliques = self._get_larger_than_stabilizer_cocliques()
        self.minimally_transitive = self._get_is_min_trans()
        self.minimally_transitive_subgroups = self._get_minimally_transitive_subgroups()
        self.is_abelian = group.is_abelian()
        self.is_nilpotent = group.is_nilpotent()
        self.is_primitive = group.is_primitive()

    def _get_derangement_classes(self):
        derangement_classes = []
        for conjugacy_class in self.conjugacy_classes:
            representative = conjugacy_class[0]
            if not Permutation(representative).fixed_points():
                derangement_classes.append((representative, len(list(conjugacy_class))))
        return derangement_classes

    def _get_eigenvalues(self):
        characters = self.characters
        derangement_classes = self.derangement_classes
        identity = self.identity
        def computation(character):
            eigenvalue_sum = 0
            eigenvalue_factor = (1/character(identity))
            for derangement_class in derangement_classes:
                character_value = character(derangement_class[0])
                eigenvalue_sum += derangement_class[1] * character_value
                # derangement_class[1] is size of its conj. class
            eigenvalue = eigenvalue_factor * eigenvalue_sum
            return (int(eigenvalue), int((character(identity) ** 2)))
        results = list(map(computation, characters))
        # loop through and consolidate multiplicities
        eigenvalues_with_multiplicities = [(k, sum(map(itemgetter(1), g)))
          for k, g in groupby(sorted(results, key=itemgetter(0), reverse=True), key=itemgetter(0))]

        return eigenvalues_with_multiplicities

    def _get_eigenvalues_subgroup(self, subgroup):
        characters = subgroup.irreducible_characters()
        conjugacy_classes = subgroup.conjugacy_classes()
        identity = subgroup.one()

        derangement_classes = []
        for conjugacy_class in conjugacy_classes:
            representative = conjugacy_class[0]
            if not Permutation(representative).fixed_points():
                derangement_classes.append(conjugacy_class)
        

        eigenvalues = []
        for character in characters:
            eigenvalue_sum = 0
            eigenvalue_factor = (1/character(identity))
            for derangement_class in derangement_classes:
                representative = derangement_class[0]
                character_value = character(representative)
                eigenvalue_sum += len(derangement_class) * character_value

            eigenvalue = eigenvalue_factor * eigenvalue_sum
            eigenvalues.append(eigenvalue)

        return eigenvalues

    def _get_max_eigenvalue(self):
        return self.eigenvalues[0][0]

    def _get_min_eigenvalue(self):
        return self.eigenvalues[-1][0]
    
    def _get_n_cliques(self):                                    
        subgroups = self.subgroups
        n_subgroups = [subgroup for subgroup in subgroups if subgroup.order() == self.degree]

        n_cliques = []
        for subgroup in n_subgroups:
            eigenvalues = self._get_eigenvalues_subgroup(subgroup)
            minimum_eigenvalue = gap.Minimum(eigenvalues)

            if minimum_eigenvalue == -1:
                n_cliques.append(subgroup)

        return n_cliques


    def _get_stabilizer_sized_cocliques(self):
        subgroups = self.subgroups
        stabilizer_sized_subgroups = [subgroup for subgroup in subgroups if subgroup.order() == self.size_of_stabilizer]

        stabilizer_sized_cocliques = []
        for subgroup in stabilizer_sized_subgroups:
            eigenvalues = self._get_eigenvalues_subgroup(subgroup)
            eigenvalue_is_zero = [eigenvalue == 0 for eigenvalue in eigenvalues]

            if all(eigenvalue_is_zero):
                stabilizer_sized_cocliques.append(subgroup)
        
        return stabilizer_sized_cocliques

    def _get_larger_than_stabilizer_cocliques(self):
        subgroups = self.subgroups
        larger_than_stabilizer_subgroups = [subgroup for subgroup in subgroups if subgroup.order() > self.size_of_stabilizer]

        larger_than_stabilizer_cocliques = []
        for subgroup in larger_than_stabilizer_subgroups:
            eigenvalues = self._get_eigenvalues_subgroup(subgroup)
            eigenvalue_is_zero = [eigenvalue == 0 for eigenvalue in eigenvalues]

            if all(eigenvalue_is_zero):
                larger_than_stabilizer_cocliques.append(subgroup)
        
        return larger_than_stabilizer_cocliques

        # size of largest gives lower bound on int. density

    def _get_is_min_trans(self):
        if self.transitivity != 0:
            indices = gap.MinimalTransitiveIndices(self.degree)
            number = gap.TransitiveIdentification(self.group)
            return (number in indices)
        return False

    def _get_minimally_transitive_subgroups(self):
        if (self.minimally_transitive):
            return []
        else:
            x = gap.MinimalTransitiveIndices(int(self.degree))
            minimal_transitive_subgroups = []

            for subgroup in self.subgroups:
                orbit_list = gap.Orbit(subgroup, 1)
                desired = list(range(1, self.degree + 1))

                if (desired == sorted(orbit_list)):
                    subgroup_index = gap.TransitiveIdentification(subgroup) 
                    #this line may take a severely long time to compute
                    #and we are not entirely sure how long yet.

                    if (subgroup_index in x):
                        minimal_transitive_subgroups.append(subgroup_index)

            unique_min_trans = set(minimal_transitive_subgroups)
            minimal_transitive_subgroups = sorted(list(unique_min_trans))

        return minimal_transitive_subgroups
