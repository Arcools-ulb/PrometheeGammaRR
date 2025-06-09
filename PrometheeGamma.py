from DataTabModels import DataTabModel
from ResultTabModel import *


class PrometheeGamma:
    """
    This class is the main class of the PROMETHEE Gamma method. It computes different matrix needed,
    such as the gamma matrix and the matrix of results

    Attributes
    ----------
    matrixGamma : list of list of float
        the gamma matrix of the PROMETHEE Gamma method
    matrixI : list
        the indifference matrix
    matrixJ : list
        the incompatibility matrix
    matrixP : list
        the preference matrix
    matrixResults : list
        the result matrix of the method
    dataTabModel : DataTabModel
        the model that contains the data of the method
    resultTabModel : ResultTabModel
        the model that contain the parameters needed for the method
    """

    def __init__(self) -> None:
        self.matrixGamma = []
        self.matrixI = []
        self.matrixJ = []
        self.matrixP = []
        self.matrixResults = []
        self.matrixGammaWeight = []
        self.dataTabModel = None
        self.resultTabModel = None


    def setDataTabModel(self, model:DataTabModel) -> None:
        """Set the dataTabModel

        Parameters
        ----------
        model : DataTabModel
            the new dataTabModel
        """
        self.dataTabModel = model


    def setResultTabModel(self, model:ResultTabModel) -> None:
        """Set the resultTabModel

        Parameters
        ----------
        model : ResultTabModel
            the new resultTabModel
        """
        self.resultTabModel = model


    def isComputed(self) -> bool:
        """Test if the results was already computed or not

        Returns
        -------
        bool
            True if the results was computed, False if not
        """
        return len(self.matrixGamma) > 0


    def getMatrixGamma(self) -> list:
        """Return the Gamma Matrix

        matrixGamma = [ [γ11, γ12, ... , γ1j, ... , γ1n] \n
                         [γ21, γ22, ... , γ2j, ... , γ2n] \n
                         [...] \n
                         [γi1, γi2, ... , γij , ..., γin] \n
                         [...] \n
                         [γn1, γn2, ... , γnj, ... , γnn] ] \n
        where n = number of alternatives ; i and j the index of the matrix

        Returns
        -------
        list of list of float
            the Gamma matrix
        """
        return self.matrixGamma

    def getMatrixGammaWeight(self) -> list:
        """Return the Gamma Matrix

        matrixGamma = [ [γ11, γ12, ... , γ1j, ... , γ1n] \n
                         [γ21, γ22, ... , γ2j, ... , γ2n] \n
                         [...] \n
                         [γi1, γi2, ... , γij , ..., γin] \n
                         [...] \n
                         [γn1, γn2, ... , γnj, ... , γnn] ] \n
        where n = number of alternatives ; i and j the index of the matrix

        Returns
        -------
        list of list of float
            the Gamma matrix
        """
        return self.matrixGammaWeight
    

    def getMatrixResults(self) -> list:
        """Return the matrix of results

        example with 3 alternatives, a1, a2 and a3: \n
        matrixResults = [ [a1 I a1, a1 P a2, a1 I a3] \n
                         [a1 P a2, a2 I a2, a2 J a3] \n
                         [a3 I a1, a3 J a2, a3 I a3] ]

        Returns
        -------
        list of list of str
            the result matrix of PROMETHEE Gamma method
        """
        return self.matrixResults

    
    def clearMatrixIJP(self) -> None:
        """Clear the matrixI, matrixJ and matrixP
        """
        self.matrixI.clear()
        self.matrixJ.clear()
        self.matrixP.clear()


    def computeGammaMatrix(self) -> None:
        """Compute the gamma matrix
        """
        self.matrixGamma.clear()
        self.dataTabModel.computeCriterionDependentValues()
        nbAlternatives = self.dataTabModel.getNumberOfAlternatives()
        nbCriteria = self.dataTabModel.getNumberOfCriteria()
        for i in range(nbAlternatives):
            self.matrixGamma.append([])
            self.matrixGammaWeight.append([])
            for j in range(nbAlternatives):
                gamma_ij = 0.0
                weight = []
                for k in range(nbCriteria):
                    tmp_gamma_ij, tmp_weight = self.dataTabModel.getGamma_ij_Criteria_k(i=i, j=j, criterion=k)
                    gamma_ij += tmp_gamma_ij
                    weight.append(tmp_weight)
                self.matrixGamma[i].append(gamma_ij)
                self.matrixGammaWeight[i].append(weight)


    def computeMatrixI(self) -> None:
        """Compute the matrixI
        """
        self.matrixI.clear()
        Ti = self.resultTabModel.getTi_float()
        for i in range(len(self.matrixGamma)):
            self.matrixI.append([])
            for j in range(len(self.matrixGamma)):
                gij = self.matrixGamma[i][j]
                gji = self.matrixGamma[j][i]
                self.matrixI[i].append(Ti - max(gij, gji))


    def computeMatrixJ(self) -> None:
        """Compute the matrixJ
        """
        self.matrixJ.clear()
        Tj = self.resultTabModel.getTj_float()
        for i in range(len(self.matrixGamma)):
            self.matrixJ.append([])
            for j in range(len(self.matrixGamma)):
                gij = self.matrixGamma[i][j]
                gji = self.matrixGamma[j][i]
                self.matrixJ[i].append(min(gij, gji) - Tj)


    def computeMatrixP(self) -> None:
        """Compute the matrixP
        """
        self.matrixP.clear()
        Pf = self.resultTabModel.getPf_float()
        for i in range(len(self.matrixGamma)):
            self.matrixP.append([])
            for j in range(len(self.matrixGamma)):
                gij = self.matrixGamma[i][j]
                gji = self.matrixGamma[j][i]
                if Pf == 0:
                    self.matrixP[i].append(0)
                else:
                    self.matrixP[i].append((gij - gji)/Pf)


    def ComputeMatrixIJP(self) -> None:
        """Compute matrixI, matrixJ and matrixP
        """
        self.computeMatrixI()
        self.computeMatrixJ()
        self.computeMatrixP()


    def computeMatrixResults(self) -> None:
        """Compute the result matrix

        ai Iγ aj ⇔ Iij ≥ max(Pij , Pji, Iij , Jij ), \n
        ai Jγ aj ⇔ Jij ≥ max(Pij , Pji, Iij , Jij ), \n
        ai Pγ aj ⇔ Pij ≥ max(Pij , Pji, Iij , Jij ), \n
        aj Pγ ai ⇔ Pji ≥ max(Pij , Pji, Iij , Jij )
        """
        self.matrixResults.clear()
        alternativesName = self.dataTabModel.getAlternativesName()
        self.resultTabModel.initScores(alternativesName)
        nbAlternatives = len(alternativesName)
        for i in range(nbAlternatives):
            self.matrixResults.append([])
            for j in range(nbAlternatives):
                Iij = self.matrixI[i][j]
                Jij = self.matrixJ[i][j]
                Pij = self.matrixP[i][j]
                Pji = self.matrixP[j][i]
                ai = alternativesName[i]
                aj = alternativesName[j]
                if(Iij >= max(Pij , Pji, Iij , Jij )):
                    self.matrixResults[i].append(ai + " I " + aj)
                elif(Jij >= max(Pij , Pji, Iij , Jij )):
                    self.matrixResults[i].append(ai + " J " + aj)
                elif(Pij >= max(Pij , Pji, Iij , Jij )):
                    self.matrixResults[i].append(ai + " P " + aj)
                    self.resultTabModel.incrementScore(ai)
                elif(Pji >= max(Pij , Pji, Iij , Jij )):
                    self.matrixResults[i].append(aj + " P " + ai)


    def computeAll(self) -> None:
        """Compute all matrix
        """
        self.computeGammaMatrix()
        self.ComputeMatrixIJP()
        self.computeMatrixResults()


    def reset(self):
        """Reset the model, i.e. clear all matrixes
        """
        self.clearMatrixIJP()
        self.matrixGamma.clear()
        self.matrixResults.clear()