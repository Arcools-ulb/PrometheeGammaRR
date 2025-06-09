from customtkinter import IntVar
import random
from sys import maxsize

from .Range import Range, RangeI, RangeJ
from .Search import ItSearch
from ..PrometheeGamma import PrometheeGamma


INFINITY = maxsize
"""The infinity (in computer terms)"""

class PreferenceLearning:
    """
    A class to implement a preference learning algorithm

    This algorithm is based on a small number of pairwise comparisons submitted to the user. 
    Based on the user's responses, an iterative search determines the possible values for the 3 parameters specific to the PROMETHEE Gamma method.
    
    Attributes
    ----------
    prometheeGamma : PrometheeGamma
        the PROMETHEE Gamma model
    master : CTkFrame
        the master frame
    alternatives : list of Alternative
        the list of all aternatives
    questions : list
        a list of peer comparisons: [(a1, a2, preference), (a1, a2, preference), ... ]
    listOfPairs : list of tuple of Aternative
        a list of pair of alternatives: [(a1, a2), (a1, a2), ... ]
    search : Search
        the instance of class for search method
    """


    def __init__(self, master, prometheeGamma:PrometheeGamma) -> None:
        self.prometheeGamma = prometheeGamma
        self.master = master
        self.alternatives = []
        self.questions = []
        self.listOfPairs = []
        self.search = ItSearch()

    
    def setAlternatives(self, alternatives:list):
        """Set the alternatives list

        Parameters
        ----------
        alternatives : list of Alternative
            a list of alternatives
        """
        self.alternatives = alternatives


    def selectFirstQuestion(self):
        """Select the first peer comparison at random

        Returns
        -------
        tuple of Alternative and IntVar
            a peer comparison (the 2 alternatives, and the value for relation between the 2 alternatives)
        """
        self.questions.clear()
        self.listOfPairs.clear()
        while True:
            a1 = random.choice(self.alternatives)
            a2 = random.choice(self.alternatives)
            if a1 != a2:
                self.listOfPairs.append((a1, a2))
                break
        question = (a1, a2, IntVar(master=self.master, value=5))
        self.questions.append(question)
        self.itSearch(True)
        return question
    

    def selectNextQuestion(self):
        """Select a non already selected peer comparison at random

        Returns
        -------
        tuple of Alternative and IntVar
            a peer comparison (the 2 alternatives, and the value for relation between the 2 alternatives)
        """
        while True:
            a1 = random.choice(self.alternatives)
            a2 = random.choice(self.alternatives)
            if a1 != a2 and ((a1, a2) not in self.listOfPairs) and ((a2, a1) not in self.listOfPairs):
                self.listOfPairs.append((a1, a2))
                break
        question = (a1, a2, IntVar(master=self.master, value=5))
        self.questions.append(question)
        self.itSearch(False)
        return question
    

    def itSearch(self, rst:bool) -> None:
        """Make an iteration of iterative search or or restart the whole search if rst is True

        Parameters
        ----------
        rst : bool
            if True, restart the search, otherwise make an iteration of iterative search
        """
        if rst:
            self.search = ItSearch()
            for q in self.questions:
                (rI, rJ, p) = self.computeRangeOfThresholdsForOneQuestion(q)
                self.search.addPair(rI, rJ, p)
                self.search.update()
        else:
            (rI, rJ, p) = self.computeRangeOfThresholdsForOneQuestion(self.questions[-1])
            self.search.addPair(rI, rJ, p)
            self.search.update()
    

    def getResults(self):
        """Return the results of the search

        Returns
        -------
        tuple of float
            the results of the search, i.e. (Imin, Imax, Jmin, Jmax, Pmin, Pmax)
        """
        return self.search.getState()
    

    def computeRangeOfThresholdsForOneQuestion(self, question:tuple):
        """Compute a range of possible thresholds for one peer comparison

        Parameters
        ----------
        question : tuple
            a peer comparison
        """
        (a1_t, a2_t, pref) = question
        a1 = self.alternatives.index(a1_t)
        a2 = self.alternatives.index(a2_t)
        preference = pref.get()
        if preference == 2:
            preference = 1 # No difference for parameters T_I, T_J and P_f
        matrixGamma = self.prometheeGamma.getMatrixGamma()
        Pmin = 1
        Pmax = INFINITY
        rI = None
        rJ = None
        if preference == 0: # indifference
            rI = RangeI(matrixGamma[a1][a2], matrixGamma[a2][a1], Pmax, Pmin) # T_I >= [Imin, Imax]
            rJ = rI
        elif preference == 1: # preference
            rI = RangeI(matrixGamma[a1][a2], matrixGamma[a2][a1], Pmax, Pmin) # T_I <= [Imin, Imax]
            rJ = RangeJ(matrixGamma[a1][a2], matrixGamma[a2][a1], Pmax, Pmin) # T_J >= [Jmin, Jmax]
        elif preference == -1: # incomparability
            rJ = RangeJ(matrixGamma[a1][a2], matrixGamma[a2][a1], Pmax, Pmin) # T_J <= [Jmin, Jmax]
            rI = rJ
        return (rI, rJ, preference)