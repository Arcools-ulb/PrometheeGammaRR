from sys import maxsize
from math import sqrt

from ..Range import RangeI, RangeJ, Range


INFINITY = maxsize
"""The infinity (in computer terms)"""


class ItSearch:
    """
    A class to make search and determine the value of PROMETHEE Gamma method parameters

    Attributes
    ----------
    Imin : float
        lowest possible value for I
    Imax : float
        highest possible value for I
    Jmin : float
        lowest possible value for J
    Jmax : float
        highest possible value for J
    Pmin : float
        lowest possible value for P
    Pmax : float
        highest possible value for P
    listOfRangeI : list of Range
        the list of range of I computed from user answers
    listOfRangeJ : list of Range
        the list of range of J computed from user answers
    listOfPreference : list of int
        the list of preference indicators from user answers
    """

    def __init__(self, Imin=0.0, Imax=1.0, Jmin=0.0, Jmax=1.0, Pmin=1.0, Pmax=INFINITY) -> None:
        """
        Parameters
        ----------
        Imin : float, optional
            lowest possible value for T_I (default is 0.0)
        Imax : float, optional
            highest possible value for T_I (default is 1.0)
        Jmin : float, optional
            lowest possible value for T_J (default is 0.0)
        Jmax : float, optional
            highest possible value for T_J (default is 1.0)
        Pmin : float, optional
            lowest possible value for P_f (default is 1.0)
        Pmax : float, optional
            highest possible value for P_f (default is INFINITY)
        """
        self.Imin = Imin
        self.Imax = Imax
        self.Jmin = Jmin
        self.Jmax = Jmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.listOfRangeI = []
        self.listOfRangeJ = []
        self.listOfPreference = []


    def getState(self) -> tuple:
        """Return the current state of the search

        Returns
        -------
        tuple of float
            the results of the search, i.e. (Imin, Imax, Jmin, Jmax, Pmin, Pmax)
        """
        return (self.Imin, self.Imax, self.Jmin, self.Jmax, self.Pmin, self.Pmax)


    def addPair(self, rI, rJ, preference:int):
        """Add a pair of I range and J range and the related preference indicator

        Parameters
        ----------
        rI : RangeI
            a range of I values
        rJ : RangeJ
            a range of J values
        preference : int
            a preference indicator (1 for preference, 0 for indifference and -1 for incomparability)
        """
        self.listOfRangeI.append(rI)
        self.listOfRangeJ.append(rJ)
        self.listOfPreference.append(preference)


    def update(self):
        """Update the search state
        """
        pref = self.listOfPreference[-1]
        rI = self.listOfRangeI[-1]
        rJ = self.listOfRangeJ[-1]
        if pref == 0:
            self.verifyPmin(rI)
            # T_I >= [Imin, Imax] ; T_J >= T_I
            self.resolveIndifference(rI)
        elif pref == -1:
            self.verifyPmin(rJ)
            # T_J <= [Jmin, Jmax] ; T_I <= T_J
            self.resolveIncomparability(rJ)
        elif pref == 1:
            # T_I <= [Imin, Imax] ; T_J >= [Jmin, Jmax] ; T_J >= T_I
            self.resolvePreference(rI, rJ)
    

    def verifyPmin(self, r):
        """Set Pmin if neccessary for the new range

        Parameters
        ----------
        r : RangeI or RangeJ
            the range of the new pair comparison
        """
        pmin = r.getPmin()
        if pmin > self.Pmin and pmin <= self.Pmax:
            self.Pmin = pmin
        elif pmin > self.Pmax:
            self.resolveConflict()


    def resolveIndifference(self, rI:RangeI):
        """Resolve indifference

        T_I >= [i_min, i_max] ; T_J >= T_I

        we have a known [Imin, Imax] from previous pair(s) or strarting value ([0,1] in this case)

        And we have rI = [i_min, i_max] from the new pair. i_min depends of Pmax and i_max depends of Pmin

        [Imin, Imax] >= [i_min, i_max] <=> Imin >= i_max

        we have to compare the 2 range and resolve any conflict

        In case of conflict, we look for the values of the 3 parameters that will minimize the mean conflict in the least squares sense.

        Parameters
        ----------
        rI : RangeI
            a range of I values
        """
        i_min = rI.getValForP(self.Pmax)
        i_max = rI.getValForP(self.Pmin)
        
        if self.Imin >= i_max:
            # There is nothing to do
            return
        elif self.Imax < i_min: # incompatible
            self.resolveConflict()
            return
        elif self.Imax >= i_max:
            if self.Jmax >= i_max:
                self.Imin = i_max # update Imin
                self.Jmin = max(self.Imin, self.Jmin)
            elif self.Jmax >= i_min:
                self.Imin = self.Jmax # update Imin
                self.Jmin = self.Jmax
                # Increase Pmin to decrease i_max such that i_max = Imin
                self.increasePminI(rI)
            else:
                self.resolveConflict()
                return
        elif self.Imax >= i_min and self.Imin >= i_min:
            # Increase Pmin to decrease i_max such that i_max = Imin
            self.increasePminI(rI)
        elif self.Imax >= i_min and self.Imin < i_min:
            dist = self.Imax - i_min
            if self.Jmax >= i_min + dist/2:
                self.Imin = i_min + dist/2 # update Imin
                self.Jmin = max(self.Imin, self.Jmin)
                # Increase Pmin to decrease i_max such that i_max = Imin
                self.increasePminI(rI)
            elif self.Jmax >= i_min:
                self.Imin = self.Jmax # update Imin
                self.Jmin = self.Jmax
                # Increase Pmin to decrease i_max such that i_max = Imin
                self.increasePminI(rI)
            else:
                self.resolveConflict()
                return
        
        self.resolveConstraints(True) # We check that all constraints are still satisfied, and if not, we solve the constraints


    def increasePminI(self, rI:RangeI):
        """Increase the Pmin value if possible (with I range) (if not, call resolveConflict())
        
        rI.max = x + y/Pmin

        <=> Pmin = y/(rI.max - x)

        rI.max must be equal to Imin

        Parameters
        ----------
        rI : RangeI
            a range of I values
        """
        x = rI.getX()
        y = rI.getY()
        if (self.Imin - x) <= 0:
            pmin = INFINITY
        else:
            pmin = y/(self.Imin - x)
        if pmin <= self.Pmax:
            self.Pmin = pmin
        else:
            self.resolveConflict()


    def resolveIncomparability(self, rJ:RangeJ):
        """Resolve incomparability

        T_J <= [j_min, j_max] ; I <= J

        we have a known [Jmin, Jmax] from previous pair(s) or strarting value ([0,1] in this case)

        And we have rJ = [j_min, j_max] from the new pair. j_min depends of Pmin and j_max depends of Pmax

        [Jmin, Jmax] <= [j_min, j_max] <=> Jmax <= j_min

        we have to compare the 2 range and resolve any conflict

        In case of conflict, we look for the values of the 3 parameters that will minimize the mean conflict in the least squares sense.

        Parameters
        ----------
        rJ : RangeJ
            a range of J values
        """
        j_min = rJ.getValForP(self.Pmin)
        j_max = rJ.getValForP(self.Pmax)
        
        if self.Jmax <= j_min:
            # There is nothing to do
            return
        elif self.Jmin > j_max: # incompatible
            self.resolveConflict()
            return
        elif self.Jmin <= j_min:
            if self.Imin <= j_min:
                self.Jmax = j_min # update Jmax
                self.Imax = min(self.Jmax, self.Imax) # T_I <= T_J
            elif self.Imin <= j_max:
                self.Jmax = self.Imin # update Jmax
                self.Imax = self.Imin
                # Increase Pmin to increase j_min such that j_min = Jmax
                self.increasePminJ(rJ)
            else:
                self.resolveConflict()
                return
        elif self.Jmin <= j_max and self.Jmax <= j_max:
            # Increase Pmin to increase j_min such that j_min = Jmax
            self.increasePminJ(rJ)
        elif self.Jmin <= j_max and self.Jmax > j_max:
            dist = j_max - self.Jmin
            if self.Imin <= j_max - dist/2:
                self.Jmax = j_max - dist/2 # update Jmax
                self.Imax = min(self.Jmax, self.Imax) # T_I <= T_J
                # Increase Pmin to increase j_min such that j_min = Jmax
                self.increasePminJ(rJ)
            elif self.Imin <= j_max:
                self.Jmax = self.Imin # update Imin
                self.Imax = self.Imin
                # Increase Pmin to increase j_min such that j_min = Jmax
                self.increasePminJ(rJ)
            else:
                self.resolveConflict()
                return

        self.resolveConstraints(True) # We check that all constraints are still satisfied, and if not, we solve the constraints


    def increasePminJ(self, rJ:RangeJ):
        """Increase the Pmin value if possible (with J range) (if not, call resolveConflict())
        
        rJ.min = x - y/Pmin

        <=> Pmin = y/(x - rJ.min)

        rJ.min must be equal to Jmax

        Parameters
        ----------
        rJ : RangeJ
            a range of J values
        """
        x = rJ.getX()
        y = rJ.getY()
        if (x - self.Jmax) <= 0:
            pmin = INFINITY
        else:
            pmin = y/(self.Jmax - x)
        if pmin <= self.Pmax:
            self.Pmin = pmin
        else:
            self.resolveConflict()


    def resolvePreference(self, rI:RangeI, rJ:RangeJ):
        """Resolve preference

        T_I <= [i_min, i_max] ; T_J >= [j_min, j_max] ; J >= I

        There are 2 parts in this resolution.

        1. T_I <= [i_min, i_max]

        we have a known [Imin, Imax] from previous pair(s) or strarting value ([0,1] in this case)

        And we have rI = [i_min, i_max] from the new pair. i_min depends of Pmax and i_max depends of Pmin

        [Imin, Imax] <= [i_min, i_max] <=> Imax <= i_min

        we have to compare the 2 range and resolve any conflict


        2. T_J >= [j_min, j_max]

        we have a known [Jmin, Jmax] from previous pair(s) or strarting value ([0,1] in this case)

        And we have rJ = [j_min, j_max] from the new pair. j_min depends of Pmin and j_max depends of Pmax

        [Jmin, Jmax] >= [j_min, j_max] <=> Jmin >= j_max

        we have to compare the 2 range and resolve any conflict

        In case of conflict, we look for the values of the 3 parameters that will minimize the mean conflict in the least squares sense.

        Parameters
        ----------
        rI : RangeI
            a range of I values
        rJ : RangeJ
            a range of J values
        """
        self.ILErI(rI) # part 1
        self.JGErJ(rJ) # part 2
        self.resolveConstraints(True) # We check that all constraints are still satisfied, and if not, we solve the constraints


    def ILErI(self, rI:RangeI):
        """Resolve preference, part 1 : T_I <= rI

        T_I <= [i_min, i_max]

        we have a known [Imin, Imax] from previous pair(s) or strarting value ([0,1] in this case)

        And we have rI = [i_min, i_max] from the new pair. i_min depends of Pmax and i_max depends of Pmin

        [Imin, Imax] <= [i_min, i_max] <=> Imax <= i_min

        we have to compare the 2 range and resolve any conflict

        Parameters
        ----------
        rI : RangeI
            a range of I values
        """
        i_min = rI.getValForP(self.Pmax)
        i_max = rI.getValForP(self.Pmin)
        
        if self.Imax <= i_min:
            # There is nothing to do
            return
        elif self.Imin > i_max: # incompatible
            self.resolveConflict()
            return
        elif self.Imin <= i_min:
            self.Imax = i_min # update Imax
        elif self.Imin <= i_max and self.Imax <= i_max:
            # Decrease Pmax to increase i_min such that i_min = Imax
            self.decreasePmaxI(rI)
        elif self.Imin <= i_max and self.Imax > i_max:
            dist = i_max - self.Imin
            self.Imax = i_max - dist/2 # update Imax
            # Decrease Pmax to increase i_min such that i_min = Imax
            self.decreasePmaxI(rI)

        
    def JGErJ(self, rJ:RangeJ):
        """Resolve preference, part 2 : T_J >= rJ

        T_J >= [j_min, j_max]

        we have a known [Jmin, Jmax] from previous pair(s) or strarting value ([0,1] in this case)

        And we have rJ = [j_min, j_max] from the new pair. j_min depends of Pmin and j_max depends of Pmax

        [Jmin, Jmax] >= [j_min, j_max] <=> Jmin >= j_max

        we have to compare the 2 range and resolve any conflict

        Parameters
        ----------
        rJ : RangeJ
            a range of J values
        """
        j_min = rJ.getValForP(self.Pmin)
        j_max = rJ.getValForP(self.Pmax)
        
        if self.Jmin >= j_max:
            # There is nothing to do
            return
        elif self.Jmax < j_min: # incompatible
            self.resolveConflict()
            return
        elif self.Jmax >= j_max:
            self.Jmin = j_max # update Jmin
        elif self.Jmax >= j_min and self.Jmin >= j_min:
            # Decrease Pmax to decrease j_max such that j_max = Jmin
            self.decreasePmaxJ(rJ)
        elif self.Jmax >= j_min and self.Jmin < j_min:
            dist = self.Jmax - j_min
            self.Jmin = j_min + dist/2 # update Jmin
            # Decrease Pmax to decrease j_max such that j_max = Jmin
            self.decreasePmaxJ(rJ)


    def decreasePmaxI(self, rI:RangeI):
        """Decrease the Pmax value if possible (with I range) (if not, call resolveConflict())
        
        rI.min = x + y/Pmax

        <=> Pmax = y/(rI.min - x)

        rI.min must be equal to Imax

        Parameters
        ----------
        rI : RangeI
            a range of I values
        """
        x = rI.getX()
        y = rI.getY()
        if (self.Imax - x) <= 0:
            pmax = INFINITY
        else:
            pmax = y/(self.Imax - x)
        if pmax >= self.Pmin:
            self.Pmax = pmax
        else:
            self.resolveConflict()


    def decreasePmaxJ(self, rJ:RangeJ):
        """Decrease the Pmax value if possible (with J range) (if not, call resolveConflict())
        
        rJ.max = x - y/Pmax

        <=> Pmax = y/(x - rJ.max)

        rJ.max must be equal to Jmin

        Parameters
        ----------
        rJ : RangeJ
            a range of J values
        """
        x = rJ.getX()
        y = rJ.getY()
        if (x - self.Jmin) <= 0:
            pmax = INFINITY
        else:
            pmax = y/(self.Jmin - x)
        if pmax >= self.Pmin:
            self.Pmax = pmax
        else:
            self.resolveConflict()


    def resolveConstraints(self, s:bool):
        """Resolve the constraints

        Constaints: J >= I ; 0 <= J <= 1 ; 0 <= I <= 1 ; 1 <= P <= INFINITY

        <=> [Jmin, Jmax] >= [Imin, Imax]

        <=> Jmin >= Imin and Jmax >= Imax

        Parameters
        ----------
        s : bool
            indicator of whether resolveConflict() has already been performed or not (True if not already performed, False otherwise)
        """
        if s:
            self.verifyIndifferences()
            self.verifyIncomparabilities()

        parameters = (self.Imin, self.Imax, self.Jmin, self.Jmax)
        for param in parameters:
            if param < 0:
                param = 0.0
            if param > 1:
                param = 1.0
        if self.Pmin < 1:
            self.Pmin = 1.0
        if self.Pmax < 1:
            self.Pmax = 1.0

        if self.Imin > self.Imax:
            if s:
                self.resolveConflict()
                s = False
            else:
                self.Imin = (self.Imin + self.Imax)/2
                self.Imax = self.Imin
        if self.Jmin > self.Jmax:
            if s:
                self.resolveConflict()
                s = False
            else:
                self.Jmin = (self.Jmin + self.Jmax)/2
                self.Jmax = self.Jmin
        if self.Pmin > self.Pmax:
            if s:
                self.resolveConflict()
                s = False
            else:
                self.Pmin = (self.Pmin + self.Pmax)/2
                self.Pmax = self.Pmin
        
        if self.Jmin >= self.Imin and self.Jmax >= self.Imax:
            # There is nothing to do
            return
        elif self.Jmax < self.Imin:
            if s:
                self.resolveConflict()
            else:
                self.Jmax = (self.Jmax + self.Imin)/2
                self.Imin = self.Jmax
        elif self.Jmin < self.Imin and self.Jmax >= self.Imax:
            self.Jmin = self.Imin
        elif self.Jmin >= self.Imin and self.Jmax < self.Imax:
            self.Imax = self.Jmax
        elif self.Jmin < self.Imin and self.Jmax < self.Imax:
            self.Jmin = self.Imin
            self.Imax = self.Jmax


    def resolveConflict(self):
        """Resolve conflict
        
        we look for the values of the 3 parameters that will minimize the mean conflict in the least squares sense

        we use a hill climbing algorithm
        
        We are certain of finding the optimum value of P thanks to the simple aspect of the optimality function: 
        there is one and only one optimum
        """
        almostFound = False
        bestI = None
        bestJ = None
        avgFitness = INFINITY
        p = 1.0
        while p <= INFINITY:
            avgI, avgJ = self.averageIandJ(p)
            fitI, fitJ = self.fitnessIandJ(p, avgI, avgJ)
            fit = sqrt(fitI**2 + fitJ**2)
            if fit < avgFitness:
                avgFitness = fit
                bestI = avgI
                bestJ = avgJ
            elif fit > avgFitness and not almostFound:
                almostFound = True
                p -= 0.1
            elif fit > avgFitness and almostFound:
                p -= 0.01
                break

            if almostFound:
                p += 0.01
            else:
                p += 0.1

            p = round(p, 2)
            if p > 100.0:
                p = INFINITY
                avgI, avgJ = self.averageIandJ(p)
                fitI, fitJ = self.fitnessIandJ(p, avgI, avgJ)
                fit = sqrt(fitI**2 + fitJ**2)
                if fit < avgFitness:
                    avgFitness = fit
                    bestI = avgI
                    bestJ = avgJ
                elif fit > avgFitness:
                    p = 100.0
                break

        self.Imax = bestI
        self.Imin = bestI
        self.Jmin = bestJ
        self.Jmax = bestJ
        self.Pmin = p
        self.Pmax = p
        self.resolveConstraints(False)


    def averageIandJ(self, p:float):
        """Compute the average value of I and J for a P value p

        Parameters
        ----------
        p : float
            the value of the parameter P

        Returns
        -------
        tuple of float
            (avgI, avgJ), the average of I and the average of J
        """
        avgI = 0
        avgJ = 0
        for i in range(len(self.listOfRangeI)):
            valueI = self.listOfRangeI[i].getValForP(p)
            avgI += valueI
            valueJ = self.listOfRangeJ[i].getValForP(p)
            avgJ += valueJ
        avgI /= len(self.listOfRangeI)
        avgJ /= len(self.listOfRangeI)

        # !!! T_I <= T_J !!!
        if avgI > avgJ:
            avg = (avgI+avgJ)/2
            avgI = avg
            avgJ = avg
        return avgI, avgJ
    

    def fitnessIandJ(self, p:float, avgI:float, avgJ:float):
        """Compute the fitness of avgI and avgJ for a P value p

        Parameters
        ----------
        p : float
            the value of the parameter P
        avgI : float
            the average of I
        avgJ : float
            the average of J
        
        Returns
        -------
        float
            (fitI, fitJ), the fitness of avgI and the fitness of avgJ
        """
        fitI = 0
        fitJ = 0
        for i in range(len(self.listOfRangeI)):
            valueI = self.listOfRangeI[i].getValForP(p)
            valueJ = self.listOfRangeJ[i].getValForP(p)
            gammaValues = self.listOfRangeI[i].getGammaValues()
            gammaSum = gammaValues[0] + gammaValues[1]
            if self.listOfPreference[i] == 0:
                error = 0
                if valueI > avgI:
                    error = valueI - avgI
                if avgI + avgJ < gammaSum:
                    error = max(error, gammaSum - (avgI+avgJ))
                fitI += error**2
            elif self.listOfPreference[i] == -1:
                error = 0
                if valueJ < avgJ:
                    error = avgJ - valueJ
                if avgI + avgJ > gammaSum:
                    error = max(error, (avgI+avgJ) - gammaSum)
                fitJ += error**2
            elif self.listOfPreference[i] == 1:
                if valueI < avgI:
                    fitI += (avgI-valueI)**2
                if valueJ > avgJ:
                    fitJ += (avgJ-valueI)**2
        fitI = sqrt(fitI)
        fitJ = sqrt(fitJ)
        return fitI, fitJ
    

    def verifyIndifferences(self):
        """Verify that T_I + T_J >= gammaij + gammaji for each indifference
        """
        for i, pref in enumerate(self.listOfPreference):
            if pref == 0:
                rI = self.listOfRangeI[i]
                rI:RangeI
                gammaValues = rI.getGammaValues()
                gammaSum = gammaValues[0] + gammaValues[1]
                if self.Imin + self.Jmin >= gammaSum:
                    # ok, nothing to do
                    pass
                elif self.Imax + self.Jmax >= gammaSum:
                    # Increase Imin and Jmin for that Imin + Jmin = gammaij + gammaji
                    val = gammaSum - (self.Imin + self.Jmin)
                    if self.Imin + val/2 <= self.Imax and self.Jmin + val/2 <= self.Jmax:
                        self.Imin += val/2
                        self.Jmin += val/2
                    elif self.Imin + val/2 > self.Imax:
                        self.Jmin += val - (self.Imax-self.Imin)
                        self.Imin = self.Imax
                    elif self.Jmin + val/2 > self.Jmax:
                        self.Imin += val - (self.Jmax-self.Jmin)
                        self.Jmin = self.Jmax
                else:
                    self.resolveConflict()
                    break


    def verifyIncomparabilities(self):
        """Verify that T_I + T_J <= gammaij + gammaji for each incomparability
        """
        for i, pref in enumerate(self.listOfPreference):
            if pref == -1:
                rJ = self.listOfRangeJ[i]
                rJ:RangeJ
                gammaValues = rJ.getGammaValues()
                gammaSum = gammaValues[0] + gammaValues[1]
                if self.Imax + self.Jmax <= gammaSum:
                    # ok, nothing to do
                    pass
                elif self.Imin + self.Jmin <= gammaSum:
                    # Decrease Imax and Jmax for that Imax + Jmax = gammaij + gammaji
                    val = (self.Imax + self.Jmax) - gammaSum
                    if self.Imax - val/2 >= self.Imin and self.Jmax - val/2 >= self.Jmin:
                        self.Imax -= val/2
                        self.Jmax -= val/2
                    elif self.Imax - val/2 < self.Imin:
                        self.Jmax -= val - (self.Imax-self.Imin)
                        self.Imax = self.Imin
                    elif self.Jmax - val/2 < self.Jmin:
                        self.Imax -= val - (self.Jmax-self.Jmin)
                        self.Jmax = self.Jmin
                else:
                    self.resolveConflict()
                    break


