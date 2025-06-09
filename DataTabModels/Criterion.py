#from customtkinter import (StringVar, IntVar, DoubleVar)

from .PreferenceFunction import PreferenceFunction


class Criterion:
    """
    Class to represent a criterion of the PROMETHEE γ method and compute values relatives to its column

    Attributes
    ----------
    name : StringVar
        the criterion name
    weight : DoubleVar
        the criterion weight
    preferenceFunction : PreferenceFunction
        the preference function used for this criterion
    pi_c_matrix : list of list of float
        the matrix of pi_c_ij (pi_c_matrix[i][j] = πc_ij = Fc(dc(ai, aj)): "how much ai is preferred over aj on criterion c")
    phi_c_list : list of float
        the list of phi_c_i (phi_c_list[i] = φc(ai) = (1/n-1) * ∑_{j=1}^{n} (πc_ij - πc_ji))
    column : list of float
        list of alternatives evaluation for this criterion
    """

    def __init__(self, name:str, weight:float, pfType:int=None, p:float=None, q:float=None) -> None:
        """
        Parameters
        ----------
        name : StringVar
            the criterion name
        weight : DoubleVar
            the criterion weight
        pftype : IntVar, optional
            the type of preference function for this criterion (default is None)
        p : DoubleVar, optional
            the preference threshold in preference function for this criterion (default is None)
        q : DoubleVar, optional
            the indifference threshold in preference function for this criterion (default is None)
        """
        self.name = name
        self.weight = weight
        self.preferenceFunction = PreferenceFunction(type=pfType, p=p, q=q)  # init the preference function
        self.pi_c_matrix = []
        self.phi_c_list = []
        self.column = []


    ##################################
    #### Getters and Setters part ####
    ##################################

    def setName(self, newName:str) -> None:
        """Set the name of the criterion

        Parameters
        ----------
        newName : StringVar
            the new criterion name
        """
        self.name = newName


    def setWeight(self, w:float) -> None:
        """Set the weight of the criterion

        w : DoubleVar
            the new criterion weight
        """
        self.weight = w


    def getName_str(self) -> str:
        """Return the name of the criterion

        Returns
        -------
        str
            the criterion name (in str)
        """
        return self.name


    def getWeight_float(self) -> float:
        """Return the weight of the criterion

        Returns
        -------
        float
            the criterion weight (in float)
        """
        return self.weight


    def getName(self) -> str:
        """Return the criterion name in a StringVar object

        Returns
        -------
        StringVar
            the criterion name (in a StringVar)
        """
        return self.name
    

    def getWeight(self) -> float:
        """Return the criterion weight in a DoubleVar object

        Returns
        -------
        DoubleVar
            the criterion weight (in a DoubleVar)
        """
        return self.weight


    def getPf_int(self) -> int:
        """Return the type of the preference function

        Returns
        -------
        int
            the type of the preference function used for this criterion (in int)
        """
        return self.preferenceFunction.getType_int()
    
    
    def getPf(self) -> int:
        """Return the type of the preference function in a IntVar object

        Returns
        -------
        Intvar
            the type of the preference function used for this criterion (in a IntVar)
        """
        return self.preferenceFunction.getType()
    

    def getP(self) -> float:
        """Return the preference threshold of the preference function in a DoubleVar object

        Returns
        -------
        DoubleVar
            the preference threshold of the preference function used for this criterion (in a DoubleVar)
        """
        return self.preferenceFunction.getP()
    

    def getP_float(self) -> float:
        """Return the preference threshold of the preference function

        Returns
        -------
        float
            the preference threshold of the preference function used for this criterion
        """
        return self.preferenceFunction.getP_float()
    

    def getQ(self) -> float:
        """Return the indifference threshold of the preference function in a DoubleVar object

        Returns
        -------
        DoubleVar
            the indifference threshold of the preference function used for this criterion (in a DoubleVar)
        """
        return self.preferenceFunction.getQ()
    

    def getQ_float(self) -> float:
        """Return the indifference threshold of the preference function

        Returns
        -------
        float
            the indifference threshold of the preference function used for this criterion
        """
        return self.preferenceFunction.getQ_float()

    def get_pi_c_matrix(self):
        return self.pi_c_matrix

    def get_pi_c_matrix_ij(self,i,j):
        return self.pi_c_matrix[i][j]* self.weight

    def get_pi_c_matrix_ij_noweight(self,i,j):
        return self.pi_c_matrix[i][j]


############################################################################################################################

    ########################
    ### Computation part ###
    ########################

    def addEvaluation(self, val:float) -> None:
        """Add an evaluation

        Parameters
        ----------
        val : float
            the new evaluation
        """
        self.column.append(val)


    def clearColumn(self):
        """Clear all value in the column
        """
        self.column.clear()


    def build_pi_c_matrix(self) -> None:
        """Build the pi matrix and the phi list of the criterion

        πc_ij = Fc(dc(ai, aj)): "how much ai is preferred over aj on criterion c" \n
        Fc is the preference function 
        and dc(ai, aj) = fc(ai) - fc(aj ) is the difference between the evaluations of ai and aj on the criterion c
        """
        for i in range(len(self.column)):
            self.pi_c_matrix.append([])
            for j in range(len(self.column)):
                pi_c_ij = self.preferenceFunction.compute_preference(self.column[i] - self.column[j])
                self.pi_c_matrix[i].append(pi_c_ij)


    def build_phi_c_list(self) -> None:
        """Build the phi list that corresponds to the criterion

        phi list: [φc(a1), φc(a2), ..., φc(ai), ..., φc(an)] where φc(ai): the mono-criterion net flow of ai on criterion c
        and n the number of alternatives
        """
        for i in range(len(self.column)):
            phi_c_i = 0.0
            for j in range(len(self.column)):
                phi_c_i += (self.pi_c_matrix[i][j] - self.pi_c_matrix[j][i])
            phi_c_i *= 1/(len(self.column)-1)

            self.phi_c_list.append(phi_c_i)


    def get_phi_c_list(self) -> list:
        """Return the phi list of the criterion

        phi list: [φc(a1), φc(a2), ..., φc(ai), ..., φc(an)] where φc(ai): the mono-criterion net flow of ai on criterion c
        and n the number of alternatives

        Returns
        -------
        list
            the list of phi_c_i (phi_c_list[i] = φc(ai) = (1/n-1) * ∑_{j=1}^{n} (πc_ij - πc_ji))
        """
        return self.phi_c_list


    def get_gamma_c_ij(self, i:int, j:int):
        """Return the value of γc_ij

        γc_ij = wc · (φc(ai) - φc(aj))
        when fc(ai) > fc(aj) and γc_ij = 0 otherwise,
        where wc is the weight associated to the criterion

        Parameters
        ----------
        i : int
            the index of the alternative i in the column
        j : int
            the index of the alternative j in the column

        Returns
        -------
        float
            the value of γc_ij
        """
        val = 0.0
        weight = 0.0
        if(i < len(self.column) and j < len(self.column)):
            if(self.column[i]>self.column[j]):
                weight = self.weight
                val = weight*(self.phi_c_list[i] - self.phi_c_list[j])
        return val , weight