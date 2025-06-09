from customtkinter import (DoubleVar)

class ResultTabModel:
    """
    A class with the model of the result tab. It contains the three new parameters of Promethee Gamma.

    Attributes
    ----------
    indifferenceThreshold : DoubleVar
        the indifference threshold of the PROMETHEE Gamma method
    incomparabilityThreshold : DoubleVar
        the incomparability threshold of the PROMETHEE Gamma method
    preferenceParameter : DoubleVar
        the preference parameter of the PROMETHEE Gamma method
    scores : dict of str:int
        a dictionnary of scores for each alternative in order to rank alternatives
    """

    def __init__(self, master) -> None:
        """
        Parameters
        ----------
        master : CTkFrame
            a tkinter master frame to link the DoubleVar
        """
        # Parameters
        self.indifferenceThreshold = 0.35#DoubleVar(master=master, value=0.0)
        self.incomparabilityThreshold = 0.35#DoubleVar(master=master, value=0.0)
        self.preferenceParameter = 1#2147483647#DoubleVar(master=master, value=1.0)
        
        # Scores for rank
        self.scores = {}


    def getTi(self) -> float:
        """Return the indifference threshold in a DoubleVar object

        Returns
        -------
        DoubleVar
            the indifference threshold of the PROMETHEE Gamma method
        """
        return self.indifferenceThreshold
    

    def getTj(self) -> float:
        """Return the incomparability threshold in a DoubleVar object

        Returns
        -------
        DoubleVar
            the incomparability threshold of the PROMETHEE Gamma method
        """
        return self.incomparabilityThreshold
    

    def getPf(self) -> float:
        """Return the preference parameter in a DoubleVar object

        Returns
        -------
        DoubleVar
            the preference parameter of the PROMETHEE Gamma method
        """
        return self.preferenceParameter
    

    def getTi_float(self) -> float:
        """Return the indifference threshold

        Returns
        -------
        float
            the indifference threshold of the PROMETHEE Gamma method
        """
        return self.indifferenceThreshold
    

    def getTj_float(self) -> float:
        """Return the incomparability threshold

        Returns
        -------
        float
            the incomparability threshold of the PROMETHEE Gamma method
        """
        return self.incomparabilityThreshold
    

    def getPf_float(self) -> float:
        """Return the preference parameter

        Returns
        -------
        float
            the preference parameter of the PROMETHEE Gamma method
        """
        return self.preferenceParameter
    

    def setTi(self, val:float) -> None:
        """Set the value of the indifference threshold

        Parameters
        ----------
        val : float
            the new value of indifference threshold
        """
        self.indifferenceThreshold= val


    def setTj(self, val:float) -> None:
        """Set the value of the incomparability threshold

        Parameters
        ----------
        val : float
            the new value of incomparability threshold
        """
        self.incomparabilityThreshold=val


    def setPf(self, val:float) -> None:
        """Set the value of the preference parameter

        Parameters
        ----------
        val : float
            the new value of preference parameter
        """
        self.preferenceParameter=val#.set(val)


    def initScores(self, alternativesNames:list):
        """Init the scores at zero

        Parameters
        ----------
        alternativesNames : list of str
            the list of altrenatives names
        """
        self.scores.clear()
        for name in alternativesNames:
            self.scores[name] = 0


    def incrementScore(self, alternativeName:str):
        """Increment the score of the alternative given in parameter

        Parameters
        ----------
        alternativeName : str
            the name of the alternative of wich the score must be incremented
        """
        self.scores[alternativeName] += 1


    def getScores(self) -> dict:
        """Return the Scores

        Returns
        -------
        dict of str:int
            the dictionnary of scrores of each alternative
        """
        return self.scores
    

    def reset(self):
        """Reset the model
        """
        self.setTi(0.0)
        self.setTj(0.0)
        self.setPf(1.0)
        self.scores.clear()