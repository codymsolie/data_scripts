# Intersection Density.sage

# this file uses results from common and ekr_determiner to compute upper and lower bounds
# on the intersection density of a group, as well as an exact value if it is known.

import mariadb
import sys

class Intersection_Density:
    def __init__(self, G, graph, ekr_determiner):
        self.G = G
        self.graph = graph
        self.has_ekr = ekr_determiner.has_ekr
        self.max_wtd_eigenvalue = ekr_determiner.max_wtd_eigenvalue
        self.upper_bound = self._get_upper_bound()
        self.lower_bound = self._get_lower_bound()
        self.exact_value = self._get_exact_value()
        
        # if intersection density determines EKR, set accordingly
        if abs(float(self.exact_value) - 1) < 0.0001 and self.has_ekr == -1:
          ekr_determiner._set_ekr(1)

###### DRIVER FUNCTIONS

    def _get_upper_bound(self):
      if self.has_ekr == 1:
        return 1

      # set to worst case scenario and decrease as we go.
      min_ub = float(self.G.order / 2)

      if not self.G.minimally_transitive:
        try:
          conn = mariadb.connect(
            user = "int_dens",
            password = "dbpass",
            host = "localhost",
            database = "intersection_density"
          )
          print("Connected to MariaDB!\n")
        except mariadb.Error as e:
          print(f"Error connecting to MariaDB platform: {e}")
          sys.exit(1)
        cursor = conn.cursor()

        densities = []
        for id in self.G.minimally_transitive_subgroups:
          cursor.execute("select int_dens from Minimally_Transitive where "\
                         "degree=? and gap_id=?", (int(self.G.degree), id))
          densities.append(cursor.fetchone()[0])

        if densities: #not empty
          least_density = float(min(densities))
          if least_density < min_ub:
            min_ub = least_density 
        
        if min_ub == 1:
          return 1
        

      if self.graph.is_cmp:
        return abs(self.G.min_eigenvalue) / self.G.size_of_stabilizer 
      if self.graph.is_pm_join:
        return (1 - self.G.min_eigenvalue) / 2
      
      # i am doing clique coclique last because common
      # takes a literal century to compute at times. 
      # when are these times? no idea as of now. 
      
      rb_ub = self.ub_ratio_bound()
      if rb_ub < min_ub:
        min_ub = rb_ub
      if self.check_int_dens(min_ub, 1):
        return 1

      else:
        nh_ub = self.ub_no_homomorphism()
        if nh_ub < min_ub:
          min_ub = nh_ub
        if self.check_int_dens(min_ub, 1):
          return 1

        else:
          cc_ub = self.ub_clique_coclique()
          if cc_ub < min_ub:
            min_ub = cc_ub
          if self.check_int_dens(min_ub, 1):
            return 1
          
      return min_ub  

    def _get_lower_bound(self):
      if self.has_ekr == 1:
        return 1
      if self.graph.is_cmp:
        return abs(self.G.min_eigenvalue) / self.G.size_of_stabilizer 
      if self.graph.is_pm_join:
        return (1 - self.G.min_eigenvalue) / 2
      if self.upper_bound == 1:
        return 1
      return self.lb_larger_than_stabilizer_cocliques() #this is our only lower bound

    def _get_exact_value(self):
      if self.has_ekr == 1:
        return 1
      if self.graph.is_cmp:
        return abs(self.G.min_eigenvalue) / self.G.size_of_stabilizer 
      if self.graph.is_pm_join:
        return (1 - self.G.min_eigenvalue) / 2
      if self.upper_bound == self.lower_bound:
          return self.upper_bound
      if self.graph.is_join:
        H = self.subgroup_by_non_derangements()
        non_der_common = Common(H)
        non_der_graph = Graph_Properties(non_der_common)
        non_der_ekr = EKR(non_der_common)
        non_der_int_dens = Intersection_Density(non_der_common, non_der_graph, non_der_ekr)
        if non_der_int_dens.upper_bound < self.upper_bound:
          self.upper_bound = non_der_int_dens.upper_bound
        if non_der_int_dens.lower_bound > self.lower_bound:
          self.lower_bound = non_der_int_dens.lower_bound
        return non_der_int_dens.exact_value

      return -1

###### HELPER FUNCTIONS

    def lb_larger_than_stabilizer_cocliques(self):
      if len(self.G.larger_than_stabilizer_cocliques) >= 1:
        max_size = int(max([
          subgroup.order() for subgroup in self.G.larger_than_stabilizer_cocliques
        ]))
        return max_size / self.G.size_of_stabilizer
      return 1 # worst possible bound if nothing is found here
    
    def ub_clique_coclique(self):
      largest_clique_size = 2 #initializing to smallest clique size 
      for subgroup in self.G.subgroups:
        if subgroup.order() > self.G.degree:
          break
        subgroup_common = Common(subgroup)
        if subgroup_common.min_eigenvalue == -1:
          if subgroup_common.max_eigenvalue == (subgroup_common.order - 1):
            if subgroup_common.order > largest_clique_size:
              largest_clique_size = subgroup_common.order
      return (self.G.degree / largest_clique_size)

    def ub_no_homomorphism(self):
      min_int_dens = (self.G.order / 2)
      if not self.G.minimally_transitive:
        try:
          conn = mariadb.connect(
            user = "int_dens",
            password = "dbpass",
            host = "localhost",
            database = "intersection_density"
          )
          print("Connected to MariaDB!\n")
        except mariadb.Error as e:
          print(f"Error connecting to MariaDB platform: {e}")
          sys.exit(1)

        cursor = conn.cursor()

        for id in self.G.minimally_transitive_subgroups:
          cursor.execute(
          "SELECT ekr,int_dens_hi FROM Groups WHERE gap_id=? AND degree=?",
          (id, int(self.G.degree)))

          row = cursor.fetchone()

          if row[0] == 1: # if group has EKR, we are done
            return 1

          if row[1] < min_int_dens:       # group does not have EKR, use int_dens_hi as new
            min_int_dens = row[1]
                                          
      return min_int_dens

    def ub_ratio_bound(self):
      if self.max_wtd_eigenvalue:
        return self.G.degree / (1 + int(round(self.max_wtd_eigenvalue)))
      return (self.G.order / 2)

    def subgroup_by_non_derangements(self):
      non_derangements = [] # holds all non-derangement elements of G
      for c in self.G.conjugacy_classes: 
        if not Permutation(c.representative()).is_derangement():
          for element in c:
             non_derangements.append(element)
      return PermutationGroup(non_derangements)
     
    def check_int_dens(self, int_dens, compare_value):
      if int_dens - compare_value < 0.00001:
        return true
      return false
