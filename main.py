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
    """
    load a existing file, parse and create objects(DTM and DTM_PRIME)
    """
    f = open("Data/Top15_SHA2014.csv", "r")
    list_name = f.readline().split(",")[1:]
    table = []
    for i in range(15):# number of alternatives. TO CHANGE
        table.append(f.readline().split(",")[1:])
    #Get all needed informations
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

def create_random_tablor():#
    """
    Generate a random evaluation table
    """
    table = []
    table_info = []
    for i in range(N):#NB of alternatives
        table.append([])

    for m in range(M):#Criteriums
        max=random.randint(100,1000)
        table_info.append(["C" +str(m) ,max])
        for n in range(N):
            table[n].append(round(random.uniform(float(max), float(max*2)),2))
    return table, table_info



def calculeWeight():#Generate a set weigts based on the number of criteriums
    """
    Generate a set of weights based on the number of criteria
    """
    list_weights = []
    max_weight = 0
    for i in range(M):
        weight = random.randint(0, 100)
        max_weight += weight
        list_weights.append(weight)
    list_weights = np.divide(list_weights, max_weight) #normalize
    return list_weights

def setupDataTabModel(DTM,DTM_PRIME,table,table_info, weight, type_function):
    """
    Create the objects(DTM and DTM_PRIME) that contains all needed informations
    """
    for i in range(len(table)):
        DTM.createAlternative(None,"A"+ str(i), table[i])
        if i < len(table)-1:
            DTM_PRIME.createAlternative(None, "A" + str(i), table[i])
    list_name =[]
    list_weights = weight
    list_type_fonction = []
    list_crit_P=[]
    list_crit_Q=[]

    for i in range(len(table_info)): #depending on the used fonction, setup the parameters P and Q.
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




"""
Set of function to plot result.
"""

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
    plt.xlim(0, 1.25)
    plt.ylim(0, 1.25)
    plt.xlabel(r"$\gamma_{ij}$")
    plt.ylabel(r"$\gamma_{ji}$")
    # Display the plot
    plt.show()

def plot_gammaij_gamma_ji_with_direction(gamma_ij, gamma_ji ,gamma_ij_PRIME, gamma_ji_PRIME,color):

    plt.plot([0, 0.35], [0.175, 0.35], color='black', linewidth=1.5, label='Reference Line')
    plt.plot([0.175, 0.35], [0, 0.35], color='black', linewidth=1.5, label='Reference Line')
    plt.plot([0.35, 1.175], [0.35, 2], color='black', linewidth=1.5, label='Reference Line')
    plt.plot([0.35, 2], [0.35, 1.175], color='black', linewidth=1.5, label='Reference Line')
    for i in range(len(gamma_ij_PRIME)):
        plt.plot([gamma_ij[i], gamma_ij_PRIME[i]], [gamma_ji[i], gamma_ji_PRIME[i]], color=color[i], linewidth=0.5, label='Reference Line')

    # Adding labels and title
    plt.xlim(0, 1.25)
    plt.ylim(0, 1.25)
    plt.xlabel(r"$\gamma_{ij}$")
    plt.ylabel(r"$\gamma_{ji}$")
    # Display the plot
    plt.show()


def main():
    #Calculation about rank reversal is computed and if is confirmed after one alternative less.
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
    type_to_analyse = [1]#, 4, 5] #List of fonction to use
    diff_borne_gamma = [[] for _ in range(len(type_to_analyse))]
    borne = [[] for _ in range(len(type_to_analyse))]
    gamma = [[] for _ in range(len(type_to_analyse))]
    while (x != 50):
        print(x)
        table, table_info = create_random_tablor()
        weight = calculeWeight()
        for i in range(len(type_to_analyse)):
            #Classical computation to get results


            DTM, DTM_PRIME = DataTabModel(), DataTabModel()
            RTM ,RTM_PRIME = ResultTabModel(None), ResultTabModel(None)
            PG,PG_PRIME = PrometheeGamma.PrometheeGamma(), PrometheeGamma.PrometheeGamma()
            DTM, DTM_PRIME = setupDataTabModel(DTM,DTM_PRIME, table, table_info, weight,type_to_analyse[i])
            # DTM is the first table
            PG.setDataTabModel(DTM)
            PG.setResultTabModel(RTM)
            PG.computeAll()
            # DTM_PRIME is the table with one less alternative
            PG_PRIME.setDataTabModel(DTM_PRIME)
            PG_PRIME.setResultTabModel(RTM_PRIME)
            PG_PRIME.computeAll()

            # Take the result and compute the differences between 2 tables
            RR=PrometheeGamma_rankReversal(N,PG, PG_PRIME,RTM.getTi(),RTM.getTj(),RTM.getPf())
            RR.isRRCalculated()

            #get datas from the actual worktable and store such that it will be used to plot some behaviours
            data_diff_borne_deplaReel, data_deplReel, data_borne = RR.get_data_to_plot()
            diff_borne_gamma[i] += data_diff_borne_deplaReel
            gamma[i] += data_deplReel
            borne[i] += data_borne
            RRDistanceBorne.extend(RR.borneMinusRealDeplacement)
            tmp1,tmp2= RR.estimationAndCalculationOfRR()
            # If a RR has been calculed and has occured
            for k in range(len(calculedRR)):
                for l in range(len(calculedRR[k])):
                    calculedRR[k][l] += tmp1[k][l]
                    confirmedRR[k][l] += tmp2[k][l]

            # some data to have nice plot
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

            #To get the direction of the gamma during the withdrawal of an alternative
            for k in range(len(RR.PG_matrix_gamma) - 1):
                for l in range(i + 1, len(RR.PG_matrix_gamma[k]) - 1):
                    if RR.PG_matrix_gamma[k][l] > RR.PG_PRIME_matrix_gamma[k][l] and RR.PG_matrix_gamma[l][k] > RR.PG_PRIME_matrix_gamma[l][k]:
                        direction_gamma[0] += 1
                    elif RR.PG_matrix_gamma[k][l] < RR.PG_PRIME_matrix_gamma[k][l] and RR.PG_matrix_gamma[l][k] < RR.PG_PRIME_matrix_gamma[l][k]:
                        direction_gamma[1] += 1
                    else:
                        direction_gamma[2] += 1
        x += 1

    #After a number X of test, plot the information/shapes.
    print("gamma direction: ",direction_gamma)
    print_RR_stat(calculedRR,confirmedRR)
    plot_difference_gamma_estimation(diff_borne_gamma)
    label_fonction_used = ["Usual", "V-Shape", "Level"]
    for i in range(len(type_to_analyse)):
        plot_gamma_and_deplaReel(label_fonction_used[i], [gamma[i], borne[i]])
    plot_gammaij_gamma_ji(gamma_ij_PRIME, gamma_ji_PRIME,color)
    plot_gammaij_gamma_ji_with_direction(gamma_ij,gamma_ji,gamma_ij_PRIME, gamma_ji_PRIME,color)#quite long

    return 0

def main_top15(): #Same fonction as main but for existing files
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
#main_top15() #Never used



