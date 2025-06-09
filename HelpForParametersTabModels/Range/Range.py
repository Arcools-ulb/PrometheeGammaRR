import numpy

class Range:
    """
    A class to represent a range of values

    Attributes
    ----------
    min : float | Range
        the lowest value in the range or a range
    max : float | Range
        the highest value in the range or a range
    """

    def __init__(self, min=None, max=None) -> None:
        """
        Parameters
        ----------
        min : float or Range
            lowest value in the range or a range
        max : float or Range
            highest value in the range or a range
        """
        self.min = min
        self.max = max


    def getMin(self):
        """Return min, the lowest value in the range or a range
        
        Returns
        -------
        float or Range
            lowest value in the range or a range
        """
        return self.min
    

    def getMax(self):
        """Return max, the highest value in the range or a range
        
        Returns
        -------
        float or Range
            highest value in the range or a range
        """
        return self.max


    def getValMin(self):
        """Return the lowest value of the range

        if min is a Range, call this function recursively, up to obtain a float value

        Returns
        -------
        float or Range
            lowest value in the range or a range
        """
        if isinstance(self.min, Range):
            return self.min.getValMin()
        else:
            return self.min
    

    def getValMax(self):
        """Return the highest value of the range

        if max is a Range, call this function recursively, up to obtain a float value

        Returns
        -------
        float or Range
            highest value in the range or a range
        """
        if isinstance(self.max, Range):
            return self.max.getValMax()
        else:
            return self.max
    

    def setMin(self, val):
        """Set the min value

        Parameters
        ----------
        val : float or Range
            the new min value
        """
        self.min = val


    def setMax(self, val):
        """Set the max value
        
        Parameters
        ----------
        val : float or Range
            the new max value
        """
        self.max = val


    def include(self, r) -> bool:
        """Test if r is strictly in the range

        Returns
        -------
        bool
            True if r is strictly in the range, False otherwise
        """
        return self.min < r and self.max > r 
    

    def print(self):
        """Print the extreme values of the range
        """
        if isinstance(self.min, Range):
            print("Min=", self.min.getValMin())
        else:
            print("Min=", self.min)
        if isinstance(self.max, Range):
            print("Max=", self.max.getValMax())
        else:
            print("Max=", self.max)


    # operator <
    def __lt__(self, r):
        if isinstance(r, Range):
            return self.max < r.getValMin()
        else: 
            return self.max < r
        

    # operator <=
    def __le__(self, r):
        if isinstance(r, Range):
            return self.max <= r.getValMin() or (self.max <= r.getValMax() and self.min <= r.getValMin())
        else: 
            return self.max <= r or (self.max >= r and self.min <= r)
        

    # operator ==
    def __eq__(self, r):
        if isinstance(r, Range):
            return self.max == r.getValMax() and self.min == r.getValMin()
        else:
            return self.min <= r and self.max >= r
        

    # operator !=
    def __ne__(self, r):
        if isinstance(r, Range):
            return not(self.max == r.getValMax() and self.min == r.getValMin())
        else:
            return not(self.min <= r and self.max >= r)
        

    # operator >
    def __gt__(self, r):
        if isinstance(r, Range):
            return self.min > r.getValMax()
        else: 
            return self.min > r
        

    # operator >=
    def __ge__(self, r):
        if isinstance(r, Range):
            return self.min >= r.getValMax() or (self.min >= r.getValMin() and self.max >= r.getValMax())
        else: 
            return self.min >= r or (self.max >= r and self.min <= r)