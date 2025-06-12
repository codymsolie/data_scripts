# Common.sage
# 
# collects general properties of the 
# provided group
#
# called upon by Generator, with
# instances passed to EKR_Determiner
# and Intersection_Density (.sage)

class Common: 
    def __init__(self, group):
        self.group = group
        self.degree = group.degree()
        self.order = group.order()
        self.one = group.one()
        self.conjugacy_classes = group.conjugacy_classes()
        self.subgroups = group.conjugacy_classes_subgroups() #avoid double work
        self.characters = group.irreducible_characters()
        self.derangement_classes = self._get_derangement_classes()
        self.eigenvalues = self._get_eigenvalues()
        self.max_eigenvalue = self._get_max_eigenvalue()
        self.min_eigenvalue = self._get_min_eigenvalue()
        self.size_of_stabilizer = self.group.stabilizer(1).order()
        self.is_join = self._get_is_join()
        self.is_cmp = self._get_is_complete_multipartite()
        self.min_trans = self._get_min_trans()

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
          for k, g in groupby(sorted(results, key=itemgetter(0)), key=itemgetter(0))]
        return eigenvalues_with_multiplicities

    def _get_max_eigenvalue(self):
        return self.eigenvalues[-1][0]

    def _get_min_eigenvalue(self):
        return self.eigenvalues[0][0]

    def _get_is_join(self):
        return int(self.max_eigenvalue - self.min_eigenvalue) == int(self.order)

    def _get_is_complete_multipartite(self):
        if self.is_join:
            eigenvalues = set(map(itemgetter(0), self.eigenvalues)) # removes multiplicity
            desired = set([int(self.max_eigenvalue), int(0), int(self.min_eigenvalue)])
            return eigenvalues == desired
        return False
    
    def _get_is_min_trans(self):

