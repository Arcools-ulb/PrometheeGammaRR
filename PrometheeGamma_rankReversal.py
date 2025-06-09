from ResultTabModel import ResultTabModel
import math
from shapely.geometry import Point, LineString #use to calculate distance between points and segments
import numpy as np
import matplotlib.pyplot as plt
class PrometheeGamma_rankReversal:
    """
    A class with the model of the result tab. It contains the three new parameters of Promethee Gamma.

    Attributes
    ----------
    """

    def __init__(self, N, PG, PG_PRIME,Ti,Tj,Pf) -> None:
        """
        Parameters
        ----------
        PG : CTkFrame
            Receive computed PROMETHEE GAMMA
        PG_PRIME : CTkFrame
            Receive computed PROMETHEE GAMMA
        """
        self.N= N
        self.PG= PG
        self.PG_PRIME = PG_PRIME

        self.PG_matrix_gamma = PG.getMatrixGamma()
        self.PG_PRIME_matrix_gamma = PG_PRIME.getMatrixGamma()

        self.PG_RESULT = PG.getMatrixResults()
        self.PG_PRIME_RESULT = PG_PRIME.getMatrixResults()

        self.Ti = Ti
        self.Tj = Tj
        self.Pf = Pf


        self.computeBorneGamma()
        self.computeBorneGamma_forPP()
        self.isBorneRight()
        self.borneMinusRealDeplacement =[]


    def computeBorneGamma(self):
        PG_WEIGHT = self.PG_PRIME.getMatrixGammaWeight()
        self.PG_estimation = []
        self.PG_estimation_positif=[]
        self.PG_estimation_negatif = []
        for i in range(len(self.PG_PRIME_matrix_gamma)):
            estimation_positif = []
            estimation_negatif=[]
            estimation_global = []
            for j in range(len(self.PG_PRIME_matrix_gamma[i])):
                nb = 0
                #print("----------")
                #print(i, j)
                for m in range(len(PG_WEIGHT[i][j])):
                    #print(PG_WEIGHT[i][j][m])
                    #print(self.PG.dataTabModel.getAlternative(i).getEvaluation(m))
                    #print(self.PG.dataTabModel.getAlternative(j).getEvaluation(m))
                    nb += PG_WEIGHT[i][j][m] * max(1, 2 * self.PG.dataTabModel.get_pi_c_ij(i, j, m))

                nb = (nb / (self.N - 2) - self.PG_matrix_gamma[i][j]/(self.N-2))
                #print(nb)
                estimation_positif.append(nb)
                estimation_negatif.append(- self.PG_matrix_gamma[i][j]/(self.N-2))
                # Python has a problem with float rounding
                # which propagates through calculations (depending on the number of alternatives)
                # We add a small number to compensate
                estimation_global.append(max(nb+0.000000000000001, (self.PG_matrix_gamma[i][j]/(self.N-2))+0.000000000000001))
            self.PG_estimation_positif.append(estimation_positif)
            self.PG_estimation_negatif.append(estimation_negatif)
            self.PG_estimation.append(estimation_global)

    def computeBorneGamma_forPP(self):
        PG_WEIGHT = self.PG_PRIME.getMatrixGammaWeight()
        self.PG_estimation_forPP=[]

        for i in range(len(self.PG_PRIME_matrix_gamma)):
            tmp = []
            for j in range(len(self.PG_PRIME_matrix_gamma[i])):
                nb = 0
                for m in range(len(PG_WEIGHT[i][j])):
                    nb += PG_WEIGHT[i][j][m] * max(1, 2 * self.PG.dataTabModel.get_pi_c_ij(i, j, m)) #Calcule a phi_c(a_i)
                # Python has a problem with float rounding
                # which propagates through calculations (depending on the number of alternatives)
                # We add a small number to compensate
                nb = (nb / (self.N - 1)) + 0.00000000000001#((nb/(self.N - 2))-(self.PG_matrix_gamma[i][j]/(self.N - 2))) +0.0000000000000001
                tmp.append(nb)
            self.PG_estimation_forPP.append(tmp)


    def isBorneRight(self):
        result = []
        for i in range(len(self.PG_PRIME_matrix_gamma)):
            for j in range(len(self.PG_PRIME_matrix_gamma[i])):
                if -self.PG_estimation[i][j] <= self.PG_matrix_gamma[i][j] - self.PG_PRIME_matrix_gamma[i][j] <= self.PG_estimation[i][j]:
                    result.append([self.PG_matrix_gamma[i][j] - self.PG_PRIME_matrix_gamma[i][j], self.PG_estimation[i][j]])
                else:
                    print("###########################################################################################")
                    print("Gamma number: ", i, j)
                    print("Gamma", self.PG_matrix_gamma[i][j], " |GMAA PRIME: ", self.PG_PRIME_matrix_gamma[i][j])
                    print("Real Movement: ", self.PG_matrix_gamma[i][j] - self.PG_PRIME_matrix_gamma[i][j])
                    print("Old_Estimation/estimation for PP relation: ", self.PG_estimation_forPP[i][j])
                    print("New_estimation_positif: ", self.PG_estimation_positif[i][j])
                    print("New_estimation_negatif: ", self.PG_estimation_negatif[i][j])
                    print("Générale_Estimation: ", self.PG_estimation[i][j])
                    print("---------------------------------------------------------------")
                    print("Is New better than Old ?: ", self.PG_estimation[i][j] >= self.PG_estimation_forPP[i][j])
                    print("Is new good ?: ",
                          self.PG_matrix_gamma[i][j] - self.PG_PRIME_matrix_gamma[i][j] <= self.PG_estimation[i][j])
                    raise ValueError('Estimation borne is too small')
        return result

    def get_data_to_plot(self):
        data_diff_borne_deplaReel = []
        data_borne = []
        data_deplReel = []
        for i in range(len(self.PG_PRIME_matrix_gamma)):
            for j in range(len(self.PG_PRIME_matrix_gamma[i])):
                data_diff_borne_deplaReel.append(self.PG_estimation[i][j]-(abs(self.PG_matrix_gamma[i][j])- abs(self.PG_PRIME_matrix_gamma[i][j])))
                data_deplReel.append( abs(self.PG_matrix_gamma[i][j] - self.PG_PRIME_matrix_gamma[i][j]))
                data_borne.append(self.PG_estimation[i][j])

        return data_diff_borne_deplaReel,data_deplReel,data_borne



    def isRRCalculated(self):
        #loop on the relation
        for i in range(len(self.PG_PRIME_RESULT)):
            for j in range(i,len(self.PG_PRIME_RESULT[i])):
                #check if relation changed
                if self.PG_RESULT[i][j] != self.PG_PRIME_RESULT[i][j]:
                    if self.specificRankReversal(i,j) == False:#Check if RR has been predictied
                        print(self.PG_RESULT[i][j], self.PG_PRIME_RESULT[i][j])
                        print("Gammaij", self.PG_matrix_gamma[i][j],"Gammaji",self.PG_matrix_gamma[j][i])
                        print("Gamma_primeij", self.PG_PRIME_matrix_gamma[i][j],"Gamma_primeji", self.PG_PRIME_matrix_gamma[j][i])
                        print("Gamma_estimation", self.PG_estimation[i][j],"True deviation", self.PG_matrix_gamma[i][j]-self.PG_PRIME_matrix_gamma[i][j])
                        print("Gamma_prime_estimation", self.PG_estimation[j][i],"True deviation", self.PG_matrix_gamma[j][i]-self.PG_PRIME_matrix_gamma[j][i])
                        raise ValueError('A RR has not been predicted/calculated')


    def extractRelation(self,i,j):
        relation = self.PG_RESULT[i][j].split(" ")[1]
        relation_bis = self.PG_PRIME_RESULT[i][j].split(" ")[1]
        if relation == "P":
            alt1 = int(self.PG_RESULT[i][j].split(" ")[0][1:])
            alt2 = int(self.PG_RESULT[i][j].split(" ")[2][1:])
            if alt1 < alt2:
                relation = "Pi"
            else:
                relation = "Pj"
        if relation_bis == "P":
            alt1 = int(self.PG_PRIME_RESULT[i][j].split(" ")[0][1:])
            alt2 = int(self.PG_PRIME_RESULT[i][j].split(" ")[2][1:])
            if alt1 < alt2:
                relation_bis = "Pi"
            else:
                relation_bis = "Pj"
        return relation, relation_bis

    def estimationAndCalculationOfRR(self):
        calculedRR = [[0,0,0,0] for k in range(4)]  # I, J, Pi,Pj
        confirmedRR = [[0, 0, 0 , 0] for k in range(4)]  # I, J, Pi, Pj
        list=["I","J","Pi", "Pj"]
        dic = {"I":0, "J":1, "Pi":2 , "Pj":3}
        for i in range(len(self.PG_PRIME_matrix_gamma)):#estimation
            for j in range(i+1, len(self.PG_PRIME_matrix_gamma[i])):
                relation, relation_bis = self.extractRelation(i,j)
                if self.PG_PRIME_RESULT[i][j] != self.PG_RESULT[i][j]:
                    confirmedRR[dic[relation]][dic[relation_bis]] +=1
                for k in list:
                    calculedRR[dic[relation]][dic[k]]+= self.calculateRankReversaleach(i,j,relation,k)
        return calculedRR, confirmedRR


    def calculateRankReversaleach(self,i,j,relation,relation_bis):
        numberOfRR =0
        if relation == "I" and relation_bis =="J":
            numberOfRR += self.RR_IJ(i,j)
        elif relation == "I" and relation_bis =="Pi":
            numberOfRR +=self.RR_IPij(i,j)
        elif relation == "I" and relation_bis == "Pj":
            numberOfRR +=self.RR_IPji(i,j)
        elif relation == "J" and relation_bis =="I":
            numberOfRR +=self.RR_IJ(i, j)
        elif relation == "J" and relation_bis =="Pi":
            numberOfRR +=self.RR_JPij(i, j)
        elif relation == "J" and relation_bis == "Pj":
            numberOfRR +=self.RR_JPji(i, j)
        elif relation == "Pi" and relation_bis == "I":
            numberOfRR +=self.RR_IPij(i, j)
        elif relation == "Pj" and relation_bis == "I":
            numberOfRR +=self.RR_IPji(i, j)
        elif relation == "Pi" and relation_bis == "J":
            numberOfRR +=self.RR_JPij(i, j)
        elif relation == "Pj" and relation_bis == "J":
            numberOfRR +=self.RR_JPji(i, j)
        elif relation == "Pi" and relation_bis == "Pj" or \
                relation_bis == "Pi" and relation == "Pj":
            if relation == "Pi":
                numberOfRR +=self.RR_Pij_To_Pji(i,j)
            else:
                numberOfRR +=self.RR_Pji_To_Pij(i,j)
        return numberOfRR

    def specificRankReversal(self,i,j):
        relation_PG, relation_PG_PRIME =self.PG_RESULT[i][j].split(" ")[1], self.PG_PRIME_RESULT[i][j].split(" ")[1]

        if relation_PG == "P" and relation_PG_PRIME == "P":
            relation_PG, relation_PG_PRIME = self.extractRelation(i,j)
            if relation_PG == "Pi":
                return self.RR_Pij_To_Pji(i,j)
            else:
                return self.RR_Pji_To_Pij(i,j)
        elif (relation_PG == "I" and relation_PG_PRIME == "J") or \
            (relation_PG == "J" and relation_PG_PRIME == "I"):
            return self.RR_IJ(i, j)
        else:
            if relation_PG == "P":
                alt1 = int(self.PG_RESULT[i][j].split(" ")[0][1:])
                alt2 = int(self.PG_RESULT[i][j].split(" ")[2][1:])
            else:
                alt1 = int(self.PG_PRIME_RESULT[i][j].split(" ")[0][1:])
                alt2 = int(self.PG_PRIME_RESULT[i][j].split(" ")[2][1:])
            if relation_PG == "I" or relation_PG_PRIME == "I":
                if alt1 < alt2:
                    return self.RR_IPij(i, j)
                else:
                    return self.RR_IPji(i, j)
            elif relation_PG == "J" or relation_PG_PRIME == "J":
                if alt1 < alt2:
                    return self.RR_JPij(i, j)
                else:
                    return self.RR_JPji(i, j)


    def RR_PijPji_version_1(self,i,j):
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        if gammaij >gammaji:
            DeltaGamma = self.PG_estimation[i][j]
        elif gammaij == gammaji:
            DeltaGamma =max(self.PG_estimation[i][j], self.PG_estimation[j][i])
        else:
            DeltaGamma = self.PG_estimation[j][i]
        if abs(gammaij - gammaji) > abs(DeltaGamma):
            return False
        self.borneMinusRealDeplacement.append(("PiPj", abs(DeltaGamma)- abs(gammaij - gammaji)))
        return True

    def RR_Pji_To_Pij(self,i,j):
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        #borne = self.PG_estimation[i][j] + self.PG_estimation[j][i]
        if gammaij > gammaji:
            borne = self.PG_estimation_forPP[i][j]
        elif gammaij == gammaji:
            borne = max(self.PG_estimation_forPP[i][j], self.PG_estimation_forPP[j][i])
        else:
            borne = self.PG_estimation_forPP[j][i]

        gamma_point = Point(gammaij, gammaji)

        Aij = ((self.Pf * self.Ti)/(1+self.Pf),0)#tuple(np.linalg.solve(np.array([[1, -1], [-1, 1/(1+self.Pf)]]), np.array([borne, -((self.Pf*self.Ti)/(1+self.Pf))])))
        Bij = (2,((self.Pf * self.Tj) / (1 + self.Pf)) + (2 / (1 + self.Pf)))#tuple(np.linalg.solve(np.array([[1, -1], [1/(1+self.Pf),-1]]), np.array([borne, -((self.Pf*self.Tj)/(1+self.Pf))])))

        polygone_Pij = LineString([Aij,(self.Ti, self.Ti),(self.Tj, self.Tj), Bij])

        if polygone_Pij.distance(gamma_point) <= borne:
            return True
        return False


    def RR_Pij_To_Pji(self,i,j):
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        #borne = self.PG_estimation[i][j] + self.PG_estimation[j][i]

        if gammaij > gammaji:
            borne = self.PG_estimation_forPP[i][j]
        elif gammaij == gammaji:
            borne = max(self.PG_estimation_forPP[i][j], self.PG_estimation_forPP[j][i])
        else:
            borne = self.PG_estimation_forPP[j][i]



        gamma_point = Point(gammaij, gammaji)

        Aji = (0,(self.Pf * self.Ti)/(1+self.Pf))#tuple(np.linalg.solve(np.array([[1, -1], [1/(1+self.Pf),-1]]), np.array([-borne, -((self.Pf*self.Ti)/(1+self.Pf))])))
        Bji = (((self.Pf * self.Tj) / (1 + self.Pf)) + (2 / (1 + self.Pf)),2)#tuple(np.linalg.solve(np.array([[1, -1], [-1,1/(1+self.Pf)]]), np.array([-borne, -((self.Pf*self.Tj)/(1+self.Pf))])))

        polygone_Pji = LineString([Aji, (self.Ti, self.Ti), (self.Tj, self.Tj), Bji])

        if polygone_Pji.distance(gamma_point) <= borne:
            return True
        return False


    def RR_IJ(self,i,j):
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        borne = (self.PG_estimation[i][j] + self.PG_estimation[j][i])/math.sqrt((-1)**2 + (-1)**2)

        if math.sqrt((gammaij - self.Ti) ** 2 + (gammaji - self.Ti) ** 2) <= borne \
                and math.sqrt((gammaij - self.Tj) ** 2 + (gammaji - self.Tj) ** 2) <= borne:
            return True
        return False

    def RR_IJ_oldversion(self, i, j):
        gammaij=self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        borne=self.PG_estimation[i][j] + self.PG_estimation[j][i]
        if self.Tj - self.Ti > borne:
            return False
        if abs(self.Tj + self.Ti - gammaij - gammaji) > abs(borne):
            return False
        self.borneMinusRealDeplacement.append(("IJ", abs(borne) - abs(self.Tj + self.Ti - gammaij - gammaji)))
        return True





    def RR_IPij(self,i,j):
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        borne_droite = (self.PG_estimation[i][j]+ self.PG_estimation[j][i]/(self.Pf + 1))
        borne_distance = borne_droite / math.sqrt((-1) ** 2 + (1 / (self.Pf + 1)) ** 2)

        gamma_point = Point(gammaij, gammaji)
        Ax = (self.Pf * self.Ti)/(1+self.Pf)

        line = LineString([(self.Ti, self.Ti), (Ax,0)])

        if line.distance(gamma_point) < borne_distance:
            return True
        return False

    def RR_IPij_version_1(self,i,j):
        borne = (self.Pf + 1)*self.PG_estimation[i][j]+ self.PG_estimation[j][i]
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        if abs(self.Pf * self.Ti - (self.Pf+1)*gammaij+gammaji) > abs(borne):
            return False
        self.borneMinusRealDeplacement.append(("IPi", abs(borne) - abs(self.Pf * self.Ti - (self.Pf+1)*gammaij+gammaji)))
        return True

    def RR_IPij_version_2(self,i,j):
        borne_droite = (self.PG_estimation[i][j]+ self.PG_estimation[j][i]/(self.Pf + 1))
        borne_distance = (self.PG_estimation[i][j] + self.PG_estimation[j][i] / (self.Pf + 1)) / math.sqrt((-1) ** 2 + (1 / (self.Pf + 1)) ** 2)
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]

        if (abs((self.Pf * self.Ti/(self.Pf + 1)) - gammaij+ gammaji/(self.Pf+1)) <= abs(borne_droite)\
                and -gammaij/(self.Pf + 1) - gammaji + (1/(self.Pf + 1) + 1)*self.Ti >= 0)\
                or math.sqrt((gammaij - self.Ti) ** 2 + (gammaji - self.Ti) ** 2) <= borne_distance:
            self.borneMinusRealDeplacement.append(("IPi", abs(borne_droite) - abs(self.Pf * self.Ti - (self.Pf + 1) * gammaij + gammaji)))
            return True
        return False






    def RR_IPji(self,i,j):
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        borne_droite = (self.PG_estimation[i][j] / (self.Pf + 1) + self.PG_estimation[j][i])
        borne_distance = borne_droite / math.sqrt((-1) ** 2 + (1 / (self.Pf + 1)) ** 2)

        gamma_point = Point(gammaij, gammaji)

        Ay = (self.Pf * self.Ti) / (1 + self.Pf)
        line = LineString([(self.Ti, self.Ti), (0, Ay)])

        if line.distance(gamma_point) < borne_distance:
            return True
        return False

    def RR_IPji_version_2(self,i,j):
        borne_droite = (self.PG_estimation[i][j] / (self.Pf + 1) + self.PG_estimation[j][i])
        borne_distance = borne_droite / math.sqrt((-1) ** 2 + (1 / (self.Pf + 1)) ** 2)
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        if (abs((self.Pf * self.Ti / (self.Pf + 1)) + gammaij / (self.Pf + 1) - gammaji) <= abs(borne_droite) \
            and -gammaij * (self.Pf + 1) - gammaji + (self.Pf + 2) * self.Ti >= 0) \
                or math.sqrt((gammaij - self.Ti) ** 2 + (gammaji - self.Ti) ** 2) <= borne_distance:
            self.borneMinusRealDeplacement.append(
                ("IPi", abs(borne_droite) - abs(self.Pf * self.Ti - (self.Pf + 1) * gammaij + gammaji)))
            return True
        return False
    # -gammaij * (self.Pf + 1) - gammaji + (self.Pf + 2)  * self.Tj <= 0
    # (abs((self.Pf * self.Tj / (self.Pf + 1)) + gammaij / (self.Pf + 1) - gammaji ) <= abs(borne_droite)


    def RR_IPji_version_1(self,i,j):
        borne =  (self.Pf + 1) * self.PG_estimation[j][i] +self.PG_estimation[i][j]
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]

        if abs(self.Pf * self.Ti + gammaij - (self.Pf + 1) * gammaji) > abs(borne):
            return False
        self.borneMinusRealDeplacement.append(("IPj", abs(borne)- abs(self.Pf * self.Ti + gammaij - (self.Pf + 1) * gammaji)))
        return True





    def RR_JPij(self,i,j):
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        borne_droite = (self.PG_estimation[i][j] / (self.Pf + 1)) + self.PG_estimation[j][i]

        borne_distance = borne_droite / math.sqrt((-1) ** 2 + (1 / (self.Pf + 1)) ** 2)


        gamma_point = Point(gammaij, gammaji)

        By = ((self.Pf * self.Tj) / (1 + self.Pf)) + (2 / (1 + self.Pf))#2*(1+self.Pf) - self.Pf * self.Tj

        line = LineString([(self.Tj, self.Tj), (2, By)])

        if line.distance(gamma_point) < borne_distance:
            return True
        return False

    def RR_JPij_version_2(self,i,j):
        borne_droite = self.PG_estimation[i][j]/ (self.Pf + 1) + self.PG_estimation[j][i]
        borne_distance = borne_droite/ math.sqrt((-1) ** 2 + (1 / (self.Pf + 1)) ** 2)
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]

        if (abs((self.Pf * self.Tj / (self.Pf + 1)) + gammaij / (self.Pf + 1) - gammaji ) <= abs(borne_droite) \
            and -gammaij * (self.Pf + 1) - gammaji + (self.Pf + 2)  * self.Tj <= 0) \
                or math.sqrt((gammaij - self.Tj) ** 2 + (gammaji - self.Tj) ** 2) <= borne_distance:
            return True
        return False

    def RR_JPij_version_1(self,i,j):
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        borne = (self.Pf + 1) * self.PG_estimation[j][i] + self.PG_estimation[i][j]

        if abs(self.Pf * self.Tj +  gammaij - (self.Pf + 1) * gammaji) > abs(borne):
            return False
        self.borneMinusRealDeplacement.append(("JPi", abs(borne)- abs(self.Pf * self.Tj +  gammaij - (self.Pf + 1) * gammaji)))
        return True




    def RR_JPji(self,i,j):
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        borne_droite = self.PG_estimation[i][j] + self.PG_estimation[j][i] / (self.Pf + 1)
        borne_distance = borne_droite / math.sqrt((-1) ** 2 + (1 / (self.Pf + 1)) ** 2)

        gamma_point = Point(gammaij, gammaji)
        Bx = ((self.Pf * self.Tj) / (1 + self.Pf)) + (2 / (1 + self.Pf))

        line = LineString([(self.Tj, self.Tj), (Bx,2)])

        if line.distance(gamma_point) <= borne_distance:
            return True
        return False

    def RR_JPji_version_2(self,i,j):
        borne_droite = self.PG_estimation[i][j] + self.PG_estimation[j][i] / (self.Pf + 1)
        borne_distance = (self.PG_estimation[i][j] + self.PG_estimation[j][i] / (self.Pf + 1)) / math.sqrt((-1) ** 2 + (1 / (self.Pf + 1)) ** 2)
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]


        if (abs((self.Pf * self.Tj / (self.Pf + 1)) - gammaij + gammaji / (self.Pf + 1)) < abs(borne_droite) \
            and -gammaij / (self.Pf + 1) - gammaji + (1 / (self.Pf + 1) + 1) * self.Tj <= 0) \
                or math.sqrt((gammaij - self.Tj) ** 2 + (gammaji - self.Tj) ** 2) <= borne_distance:
            self.borneMinusRealDeplacement.append(("IPi", abs(borne_droite) - abs(self.Pf * self.Ti - (self.Pf + 1) * gammaij + gammaji)))
            return True

        return False


    def RR_JPji_version_1(self,i,j):
        borne = (self.Pf + 1) * self.PG_estimation[i][j]+ self.PG_estimation[j][i]
        gammaij = self.PG_matrix_gamma[i][j]
        gammaji = self.PG_matrix_gamma[j][i]
        if abs(self.Pf * self.Tj - (self.Pf + 1) * gammaij + gammaji) > abs(borne):
            return False
        self.borneMinusRealDeplacement.append(
            ("JPj", abs(borne) -abs(self.Pf * self.Tj - (self.Pf + 1) * gammaij + gammaji)))
        return True





    def reset(self):
        """Reset the model
        """
        pass