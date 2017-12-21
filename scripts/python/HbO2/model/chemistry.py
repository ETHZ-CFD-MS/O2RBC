from utilities.decorators import accept_scalar_arg


class BloodChemistry:
    """
    Equations for kinetics between hemoglobin saturation to oxygen 
    partial pressure
    """

    def __init__(self, nHill=2.64, P50=47.9):
        self.nHill = nHill
        self.P50 = P50

    @accept_scalar_arg
    def hillS(self, P):
        return pow(P.clip(min=0), self.nHill) \
             /(pow(P, self.nHill) + pow(self.P50, self.nHill))

    @accept_scalar_arg
    def hillP(self, S):
        return pow(S.clip(min=0)/(1-S), 1/self.nHill)*self.P50

    @accept_scalar_arg
    def dPdSHill(self, S):
        denominator = self.nHill*S*(1 - S)
        denominator[denominator <= 0] = 1
        return self.hillP(S)/denominator

    @accept_scalar_arg
    def dSdPHill(self, P):
        return self.nHill*pow(self.P50, self.nHill)*pow(P.clip(min=0), self.nHill-1) \
             / (pow(P.clip(min=0), self.nHill) + pow(self.P50, self.nHill))**2
