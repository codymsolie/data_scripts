# EKR.sage
#
# uses Common class data to determine
# if the EKR property holds for the group
#
# Author: Cody Antal (2021)
#
# *Common class is referred to as "G" in 
# this script.*

class EKR:
    def __init__(self, G):
        self.G = G
        self.has_ekr = -1 
        self.reasons = []

        #weightings is used (if we have some) to pass weighting information along to the EKRM determiner
        #max_wtd_eigenvalue for use in intersection density ratio bound

        self.weightings = None
        self.max_wtd_eigenvalue = None 

        if self.G.larger_than_stabilizer_cocliques:
            self.has_ekr = 0 
            self.reasons.append("A coclique of size larger than |G|/n exists")

        if self.has_ekr == -1 :
            if self._ratiobound_gives_ekr():
                self.has_ekr = 1
                self.reasons.append("Ratiobound gives EKR property")

            if self.G.n_cliques:
                self.has_ekr = 1 
                self.reasons.append("Group has a clique of size n")
            else:
                (weightings_work, weightings) = self._weighting_gives_ekr()
                if weightings_work:
                    self.has_ekr = 1 
                    self.reasons.append("A weighting on derangement classes gives EKR")

                    self.weightings = weightings
        
    
    def _ratiobound_gives_ekr(self):
        if self.G.degree - 1 == -(self.G.max_eigenvalue)/(self.G.min_eigenvalue):
            return True
        else:
            return False
    
    def _weighting_gives_ekr(self):
        coefficents = self._get_coefficents()
        open_conjugacy_classes = self._get_open_conjugacy_classes()
        
        self._create_lp(coefficents, open_conjugacy_classes)
        self._solve_lp()
        
        return self._read_sol()
    

    #helper functions

    def _get_coefficents(self):
        characters = self.G.characters
        derangement_classes = self.G.derangement_classes
        identity = self.G.identity

        coefficents = []
        for character in characters:
            coefficents_for_character = []

            coefficent_factor = 1/character(identity)
            coefficent_sum = 0
            for derangement_class in derangement_classes:
                representative = derangement_class[0]
                character_value = real(character(representative)) 
                # we ignore the the imaginary components as they will 
                # all add to 0 anyway (as long as the non-inverse
                # closed classes are weighted correctly!)

                coefficent_sum += character_value * derangement_class[1]

                coefficent = coefficent_factor * coefficent_sum
                coefficents_for_character.append(float(coefficent)) 
                # we cast to a float here to make sure the coefficent 
                # is represented as a decimal (1.41) and not as a sage Real (sqrt(2))
            
            coefficents.append(coefficents_for_character)
        
        return coefficents
    
    def _get_open_conjugacy_classes(self):
        conjugacy_classes = self.G.conjugacy_classes

        open_conjugacy_classes = []
        for conjugacy_class in conjugacy_classes:
            representative = conjugacy_class[0]
            representative_inverse = representative.inverse()

            conjugacy_class_inverse = self.G.group.conjugacy_class(representative_inverse)

            if not (conjugacy_class_inverse, conjugacy_class) in open_conjugacy_classes and conjugacy_class != conjugacy_class_inverse:
                open_conjugacy_classes.append((conjugacy_class, conjugacy_class_inverse))
        
        return open_conjugacy_classes
    


    def _create_lp(self, coefficents, open_conjugacy_classes):
        lp = open("gurobi/temp.lp", "w")
        
        lp.write("Maximize\n")
        objective_coefficents = coefficents[0]
        objective = self._create_linear_combination_string(objective_coefficents, "")
        lp.write(objective)

        lp.write("Subject To\n")
        for constraint_coefficents in coefficents[1:]:
            upper_bound = self._create_linear_combination_string(constraint_coefficents, f" <= {self.G.degree - 1}")
            lower_bound = self._create_linear_combination_string(constraint_coefficents, " >= -1")
            lp.write(upper_bound)
            lp.write(lower_bound)
        
        for (class_1, class_2) in open_conjugacy_classes:
            class_1_index = self.G.group.conjugacy_classes().index(class_1)
            class_2_index = self.G.group.conjugacy_classes().index(class_2)

            constraint = f"x{class_1_index} - x{class_2_index} = 0\n"
            lp.write(constraint)
        
        lp.write("End\n")
        lp.close()


    def _create_linear_combination_string(self, coefficents, bound):
        linear_combination = ""
        for (i, coefficent) in enumerate(coefficents):
            if coefficent < 0:
                linear_combination += f"{coefficent} x{i} "
            else:
                linear_combination += f"+{coefficent} x{i} " if i != 0 else f"{coefficent} x{i} "
            
            if i == len(coefficents) - 1:
                linear_combination += bound + "\n"
        
        return linear_combination


    def _solve_lp(self):
        os.system("gurobi_cl ResultFile=gurobi/temp.sol OutputFlag=0 LogFile=\"\" gurobi/temp.lp") 
        # this may not work on every system! (if you're trying to run any of this on windows... good luck)

    def _read_sol(self, tolerance = 0.0001):
        file = open("gurobi/temp.sol", "r")

        maximum_eigenvalue = 0
        weightings = []
        for line in file:
            if ("Objective value" in line):
                index_of_start_of_number = 20
                maximum_eigenvalue = float(line[index_of_start_of_number:])
            else:
                index_of_space = line.index(" ")
                weighting = float(line[index_of_space + 1:])
                weightings.append(weighting)

        self.max_wtd_eigenvalue = maximum_eigenvalue

        weighted_ratiobound = self.G.order/(1 + maximum_eigenvalue)
        size_of_stabilizer = self.G.size_of_stabilizer

        if abs(weighted_ratiobound - size_of_stabilizer) < tolerance:
            return (True, weightings)
        else:
            return (False, None)

    # used if intersection density determines EKR
    def _set_ekr(self, value):
      self.has_ekr = value
      return
