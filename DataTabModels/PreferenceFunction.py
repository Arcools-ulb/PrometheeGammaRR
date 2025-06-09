import numpy as np
#from customtkinter import (IntVar, DoubleVar)

class PreferenceFunction:
    """
    This class represent the preference funtion used by the Promethee Gamma method

    Attributes
    ----------
    type : IntVar
        the type of the preference function
    pc : DoubleVar
        the preference threshold in the preference function
    qc : DoubleVar
        the indifference threshold in the preference function
    """

    def __init__(self, type:int=None, p:float=None, q:float=None):
        """
        Parameters
        ----------
        type : IntVar, optional
            the type of the preference function (default is None)
        p : DoubleVar, optional
            the preference threshold in the preference function (default is None)
        q : DoubleVar, optional
            the indifference threshold in the preference function (default is None)
        """

        self.type = type
        self.pc = p
        self.qc = q


    def getP(self) -> float:
        """Return the preference threshold in a DoubleVar object

        Returns
        -------
        DoubleVar
            the preference threshold in the preference function (in a DoubleVar)
        """
        return self.pc
    

    def getP_float(self) -> float:
        """Return the preference threshold

        Returns
        -------
        float
            the preference threshold in the preference function
        """
        return self.pc
    

    def getQ(self):
        """Return the indifference threshold in a DoubleVar object

        Returns
        -------
        DoubleVar
            the indifference threshold in the preference function (in a DoubleVar)
        """
        return self.qc
    

    def getQ_float(self):
        """Return the indifference threshold

        Returns
        -------
        float
            the indifference threshold in the preference function
        """
        return self.qc
    

    def getType(self):
        """Return the type of the preference function in an IntVar object

        Returns
        -------
        IntVar
            the type of the preference function (in a IntVar)
        """
        return self.type


    def getType_int(self) -> int:
        """Return the type of the preference function

        Returns
        -------
        int
            the type of the preference function
        """
        return self.type


    def set_type(self, new_type:int, p:float=None, q:float=None) -> None:
        """Set the function

        Parameters
        ----------
        new_type : IntVar
            the new type of the preference function
        p : DoubleVar, optional
            the new preference threshold in the preference function (default is None)
        q : DoubleVar, optional
            the new indifference threshold in the preference function (default is None)
        """
        self.type = new_type
        self.pc = p
        self.qc = q


    def set_qc(self, new_qc:float) -> None:
        """Set the indifference parameter

        Parameters
        ----------
        new_qc : DoubleVar
            the new indifference threshold in the preference function
        """
        self.qc = new_qc


    def set_pc(self, new_pc:float) -> None:
        """Set the preference parameter

        Parameters
        ----------
        new_pc : DoubleVar
            the new preference threshold in the preference function
        """
        self.pc = new_pc


    def compute_preference(self, value:float) -> float:
        """Compute the preference value

        Parameters
        ----------
        value : float
            the difference between the evaluations of the alternatives i and j for the criterion that use this preference function

        Returns
        -------
        float
            the value of preference between alternatives i and j
        """
        preference = 0.0
        type = self.type
        pc = self.pc
        qc = self.qc

        if type == 1: # Usual
            if value > 0:
                preference = 1.0
        elif type == 2: # U-shape
            if value >= qc:
                preference = 1.0
        elif type == 3: # V-shape
            if value >= pc :
                preference = 1.0
            elif value > 0:
                preference = value/pc
        elif type == 4: # Level
            if value >= pc:
                preference = 1.0
            elif value >= qc:
                preference = 0.5
        elif type == 5: # Linear
            if value >= pc:
                preference = 1.0
            elif value > qc:
                preference = (value - qc)/(pc - qc)
        elif type == 6: # Gaussian
            if value > 0:
                preference = np.exp(-value**2)
        return preference