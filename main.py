import math

import random
import matplotlib.pyplot as plt
import numpy as np

import PrometheeGamma
from PrometheeGamma_rankReversal import *
from DataTabModels.DataTabModel import DataTabModel
from ResultTabModel import *

#GLOBAL VARIABLE
N = 25 #number of solutions/alternatives
M = 10   #number of criteria
X = 0

def open_top15_SHA2014(DTM,DTM_PRIME):
    f = open("Data/Top15_SHA2014.csv", "r")
    list_name = f.readline().split(",")[1:]
    table = []
    for i in range(15):# number of alternatives. TO CHANGE
        table.append(f.readline().split(",")[1:])

    list_weights = f.readline().split(",")[1:]
    list_type_fonction = f.readline().split(",")[1:]
    list_crit_P = f.readline().split(",")[1:]
    list_crit_Q = f.readline().split(",")[1:]
    f.close()
    for i in range(len(table)):
        DTM.createAlternative(None,"A"+ str(i), table[i])
        if i < len(table) -1 :
            DTM_PRIME.createAlternative(None, "A" + str(i), table[i])


    DTM.createCriteria(None, list_name, list_weights, list_type_fonction, list_crit_P, list_crit_Q)
    DTM_PRIME.createCriteria(None, list_name, list_weights, list_type_fonction, list_crit_P, list_crit_Q)
    return DTM, DTM_PRIME, table

def create_random_tablor():
    table = []
    table_info = []
    for i in range(N):
        table.append([])

    for m in range(M):
        max=random.randint(100,1000)
        table_info.append(["C" +str(m) ,max])
        for n in range(N):
            table[n].append(round(random.uniform(float(max), float(max*2)),2))
        #table[-1].append(0)
        #table[-1].append(10000)
        #table[-1].append(max + max/2)
    #print(table)
    return table, table_info



def calculeWeight():
    list_weights = []
    max_weight = 0
    for i in range(M):
        weight = random.randint(0, 100)
        max_weight += weight
        list_weights.append(weight)
    list_weights = np.divide(list_weights, max_weight)
    return list_weights #[1,0,0]

def setupDataTabModel(DTM,DTM_PRIME,table,table_info, weight, type_function):
    for i in range(len(table)):
        DTM.createAlternative(None,"A"+ str(i), table[i])
        if i < len(table)-1:
            DTM_PRIME.createAlternative(None, "A" + str(i), table[i])
    list_name =[]
    list_weights = weight
    list_type_fonction = []
    list_crit_P=[]
    list_crit_Q=[]

    for i in range(len(table_info)):
        list_name.append(table_info[i][0])
        if type_function == 4:
            list_type_fonction.append(4)
            list_crit_Q.append(1/3*table_info[i][1])
            list_crit_P.append(2/3* table_info[i][1])
        elif type_function == 5:
            list_type_fonction.append(5)
            list_crit_Q.append(0)
            list_crit_P.append(2/3*table_info[i][1])
        elif type_function == 1:
            list_type_fonction.append(1)
            list_crit_Q.append(0)
            list_crit_P.append(0)

    DTM.createCriteria(None,list_name,list_weights,list_type_fonction,list_crit_P,list_crit_Q)
    DTM_PRIME.createCriteria(None, list_name, list_weights, list_type_fonction, list_crit_P, list_crit_Q)
    return DTM, DTM_PRIME


def testBorne(PG, PG_PRIME):
    PG_matrix_gamma = PG.getMatrixGamma()
    PG_PRIME_matrix_gamma = PG_PRIME.getMatrixGamma()
    # relation
    PG_WEIGHT = PG_PRIME.getMatrixGammaWeight()
    # print(PG_WEIGHT)
    result = []
    error = 0
    for i in range(len(PG_PRIME_matrix_gamma)):
        for j in range(len(PG_PRIME_matrix_gamma[i])):
            nb = 0
            for m in range(len(PG_WEIGHT[i][j])):
                nb += PG_WEIGHT[i][j][m] * max(1, 2 * PG.dataTabModel.get_pi_c_ij(i, j, m))
            nb = nb / (N - 1)
            if -nb <= PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j] <= nb:
                result.append([PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j], nb])
            else:
                error = 1
                print(nb, PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j],
                      nb >= PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j])
                raise ValueError('Borne trop grande !!!')
    return result, error

def gammaBorneVerification(PG, PG_PRIME):
    #PG.getMatrixGamma(), PG_PRIME.getMatrixGamma(),PG.getMatrixResults(),PG_PRIME.getMatrixResults(), PG_PRIME.getMatrixGammaWeight(),PG_PRIME
    PG_matrix_gamma = PG.getMatrixGamma()
    PG_PRIME_matrix_gamma = PG_PRIME.getMatrixGamma()
    #relation
    #PG_RESULT = PG.getMatrixResults()
    #PG_PRIME_RESULT = PG_PRIME.getMatrixResults()
    PG_WEIGHT = PG_PRIME.getMatrixGammaWeight()
    #print(PG_WEIGHT)
    result = []
    error = 0
    for i in range(len(PG_PRIME_matrix_gamma)):
        for j in range(len(PG_PRIME_matrix_gamma[i])):
            if i != j:
                nb = 0
                for m in range(len(PG_WEIGHT[i][j])):
                    #print("iCIIIIC", PG.dataTabModel.get_pi_c_ij(i,j,m))
                    #print("------------------------------")
                    #print(" Poids", PG_WEIGHT[i][j])
                    #print(" Pc", PG.dataTabModel.get_pi_c_ij(i,j,m))
                    #print("max", max(1, 2 * PG.dataTabModel.get_pi_c_ij(i, j, m)))
                    nb += PG_WEIGHT[i][j][m]* max(1, 2*PG.dataTabModel.get_pi_c_ij(i,j,m))
                nb = nb/(N-1)
                print(nb,PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j] , nb >= PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j])
                if -nb <= PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j] <= nb:
                    #print(PG[i][j] - PG_PRIME[i][j], -2* PG_WEIGHT[i][j]/(N-1), -2/(N-1))
                    result.append([PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j],nb])
                else:
                    error = 1
                    print(nb,PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j] , nb >= PG_matrix_gamma[i][j] - PG_PRIME_matrix_gamma[i][j])
                    raise ValueError('Borne trop grande !!!')
    return result, error


def plotResult(result, type_to_analyse):
    number = []
    gamma_minus_gamma_prime = []
    estimation = []
    for i in range(len(result)):
        number.append(i)
        gamma_minus_gamma_prime.append(abs(result[i][0]))
        estimation.append(abs(result[i][1]))
    plt.bar(number,estimation,color='blue')
    plt.bar(number, gamma_minus_gamma_prime ,color='red')
    #plt.scatter(gamma_minus_gamma_prime, number)
    #naming
    plt.legend(['Estimation', 'Déplacement réel'])
    plt.xlabel('Déplacement de GAMMA')
    plt.show()


def plot(Value, relation):
    plt.scatter(Value, Value, color='blue')
    # plt.scatter(gamma_minus_gamma_prime, number)
    # naming
    plt.title("Relation : "+relation)
    plt.xlabel('Difference entre la borne et gamma')
    plt.ylabel('Difference entre la borne et gamma')
    plt.show()


def print_RR_stat(calculedRR, confirmedRR):
    relation = ["I", "J", "Pi", "Pj"]
    print("Calculate")
    print(relation)
    for i in range(len(calculedRR)):
        print(relation[i],calculedRR[i])
    print("effectif")
    print(relation)
    for i in range(len(confirmedRR)):
        print(relation[i], confirmedRR[i])


def plot_difference_gamma_estimation(data):
    plt.hist(data, bins=30,label=["Usual","V-Shape","Level"], edgecolor='black')

    # Adding labels and title
    plt.xlabel("Différence entre la borne et le déplacement réel")
    plt.ylabel('Fréquence')
    plt.legend()
    # Display the plot
    plt.show()

def plot_gamma_and_deplaReel(name,data):
    plt.hist(data, bins=30,label=["Réel","Borne"], edgecolor='black')

    # Adding labels and title
    plt.title("Fonction " + name)
    plt.xlabel("Déplacement")
    plt.ylabel('Fréquence')
    plt.legend()
    # Display the plot
    plt.show()

def plot_gammaij_gamma_ji(gamma_ij, gamma_ji ,color):

    plt.plot([0, 0.3], [0.15, 0.3], color='black', linewidth=1.5, label='Reference Line')
    plt.plot([0.15, 0.3], [0, 0.3], color='black', linewidth=1.5, label='Reference Line')
    plt.plot([0.3, 0.4], [0.3, 0.4], color='black', linewidth=1.5, label='Reference Line')
    plt.scatter(gamma_ij, gamma_ji ,color=color,s=1)

    # Adding labels and title
    #plt.title("Fonction ")
    plt.xlim(0, 1.25)
    plt.ylim(0, 1.25)
    plt.xlabel(r"$\gamma_{ij}$")
    plt.ylabel(r"$\gamma_{ji}$")
    #plt.legend()
    # Display the plot
    plt.show()

def plot_gammaij_gamma_ji_with_direction(gamma_ij, gamma_ji ,gamma_ij_PRIME, gamma_ji_PRIME,color):

    plt.plot([0, 0.35], [0.175, 0.35], color='black', linewidth=1.5, label='Reference Line')
    plt.plot([0.175, 0.35], [0, 0.35], color='black', linewidth=1.5, label='Reference Line')
    plt.plot([0.35, 1.175], [0.35, 2], color='black', linewidth=1.5, label='Reference Line')
    plt.plot([0.35, 2], [0.35, 1.175], color='black', linewidth=1.5, label='Reference Line')
    for i in range(len(gamma_ij_PRIME)):
        plt.plot([gamma_ij[i], gamma_ij_PRIME[i]], [gamma_ji[i], gamma_ji_PRIME[i]], color=color[i], linewidth=0.5, label='Reference Line')
    #plt.scatter(gamma_ij, gamma_ji ,color=color,s=1)

    # Adding labels and title
    #plt.title("Fonction ")
    plt.xlim(0, 1.25)
    plt.ylim(0, 1.25)
    plt.xlabel(r"$\gamma_{ij}$")
    plt.ylabel(r"$\gamma_{ji}$")
    #plt.legend()
    # Display the plot
    plt.show()


def main():
    calculedRR = [[0, 0, 0,0] for k in range(4)]  # I, J, Pi, Pj
    confirmedRR = [[0, 0, 0,0] for k in range(4)]  # I, J, Pi, Pj
    RRDistanceBorne =[]

    gamma_ij = []
    gamma_ji = []
    gamma_ij_PRIME = []
    gamma_ji_PRIME = []
    color = []

    direction_gamma = [0,0,0] #both negative, both positive, neg and pos

    x = 0
    type_to_analyse = [1]#, 4, 5]
    diff_borne_gamma = [[] for _ in range(len(type_to_analyse))]
    borne = [[] for _ in range(len(type_to_analyse))]
    gamma = [[] for _ in range(len(type_to_analyse))]
    while (x != 50):
        #print("#################NOUVELLE BOUCLE##################")
        print(x)
        table, table_info = create_random_tablor()
        weight = calculeWeight()
        for i in range(len(type_to_analyse)):
            DTM, DTM_PRIME = DataTabModel(), DataTabModel()
            RTM ,RTM_PRIME = ResultTabModel(None), ResultTabModel(None)
            PG,PG_PRIME = PrometheeGamma.PrometheeGamma(), PrometheeGamma.PrometheeGamma()
            DTM, DTM_PRIME = setupDataTabModel(DTM,DTM_PRIME, table, table_info, weight,type_to_analyse[i])
            PG.setDataTabModel(DTM)
            PG.setResultTabModel(RTM)
            PG.computeAll()


            PG_PRIME.setDataTabModel(DTM_PRIME)
            PG_PRIME.setResultTabModel(RTM_PRIME)
            PG_PRIME.computeAll()

            #plotResult(gammaBorneVerification(PG,PG_PRIME)[0], "")
            RR=PrometheeGamma_rankReversal(N,PG, PG_PRIME,RTM.getTi(),RTM.getTj(),RTM.getPf())


            RR.isRRCalculated()
            data_diff_borne_deplaReel, data_deplReel, data_borne = RR.get_data_to_plot()
            diff_borne_gamma[i] += data_diff_borne_deplaReel
            gamma[i] += data_deplReel
            borne[i] += data_borne
            #RR.loopOnRelation()
            RRDistanceBorne.extend(RR.borneMinusRealDeplacement)
            tmp1,tmp2= RR.estimationAndCalculationOfRR()

            for k in range(len(calculedRR)):
                for l in range(len(calculedRR[k])):
                    calculedRR[k][l] += tmp1[k][l]
                    confirmedRR[k][l] += tmp2[k][l]
            for k in range(len(RR.PG_matrix_gamma) - 1):
                for l in range(i+1,len(RR.PG_matrix_gamma[k]) - 1):

                    if RR.PG_RESULT[k][l] != RR.PG_PRIME_RESULT[k][l]:
                        color.append("red")
                    elif RR.specificRankReversal(k,l) == True:
                        color.append("green")
                    else:
                        color.append('blue')
                    gamma_ij.append(RR.PG_matrix_gamma[k][l])
                    gamma_ji.append(RR.PG_matrix_gamma[l][k])
                    gamma_ij_PRIME.append(RR.PG_PRIME_matrix_gamma[k][l])
                    gamma_ji_PRIME.append(RR.PG_PRIME_matrix_gamma[l][k])
            for k in range(len(RR.PG_matrix_gamma) - 1):
                for l in range(i + 1, len(RR.PG_matrix_gamma[k]) - 1):
                    if RR.PG_matrix_gamma[k][l] > RR.PG_PRIME_matrix_gamma[k][l] and RR.PG_matrix_gamma[l][k] > RR.PG_PRIME_matrix_gamma[l][k]:
                        direction_gamma[0] += 1
                    elif RR.PG_matrix_gamma[k][l] < RR.PG_PRIME_matrix_gamma[k][l] and RR.PG_matrix_gamma[l][k] < RR.PG_PRIME_matrix_gamma[l][k]:
                        direction_gamma[1] += 1
                    else:
                        direction_gamma[2] += 1
            #result, error = gammaBorneVerification(PG, PG_PRIME)
            #RR.isRelationChange()
            #print('ici',PG_PRIME.dataTabModel.get_pi_ij(0,2)/N-1)
            #print('ici', PG.dataTabModel.get_pi_ij(0, 2) / N - 1)
            #print(PG.getMatrixGamma())
        #print(PG_PRIME.getMatrixGamma())
        x += 1
    print("gamma direction: ",direction_gamma)
    print_RR_stat(calculedRR,confirmedRR)
    plot_difference_gamma_estimation(diff_borne_gamma)
    label_fonction_used = ["Usual", "V-Shape", "Level"]
    for i in range(len(type_to_analyse)):
        plot_gamma_and_deplaReel(label_fonction_used[i], [gamma[i], borne[i]])
    plot_gammaij_gamma_ji(gamma_ij_PRIME, gamma_ji_PRIME,color)
    plot_gammaij_gamma_ji_with_direction(gamma_ij,gamma_ji,gamma_ij_PRIME, gamma_ji_PRIME,color)#quite long

    return 0

def main_top15():
    calculedRR = [[0, 0, 0, 0] for k in range(4)]  # I, J, Pi, Pj
    confirmedRR = [[0, 0, 0, 0] for k in range(4)]  # I, J, Pi, Pj
    RRDistanceBorne = []


    DTM, DTM_PRIME = DataTabModel(), DataTabModel()
    DTM , DTM_PRIME, table =open_top15_SHA2014(DTM,DTM_PRIME)
    RTM, RTM_PRIME = ResultTabModel(None), ResultTabModel(None)
    PG, PG_PRIME = PrometheeGamma.PrometheeGamma(), PrometheeGamma.PrometheeGamma()
    PG.setDataTabModel(DTM)
    PG.setResultTabModel(RTM)
    PG.computeAll()

    PG_PRIME.setDataTabModel(DTM_PRIME)
    PG_PRIME.setResultTabModel(RTM_PRIME)
    PG_PRIME.computeAll()

    # plotResult(gammaBorneVerification(PG,PG_PRIME)[0], "")
    RR = PrometheeGamma_rankReversal(N, PG, PG_PRIME, RTM.getTi(), RTM.getTj(), RTM.getPf())

    RR.isRRCalculated()
    RRDistanceBorne.extend(RR.borneMinusRealDeplacement)
    tmp1, tmp2 = RR.estimationAndCalculationOfRR()

    for k in range(len(calculedRR)):
        for l in range(len(calculedRR[k])):
            calculedRR[k][l] += tmp1[k][l]
            confirmedRR[k][l] += tmp2[k][l]
    print_RR_stat(calculedRR, confirmedRR)



random.seed(1)
main()
#main_top15()



