class Func:
    def __init__(self, f) -> None:
        self.func = f
    
    def eval(self, x) -> float:
        return self.func(x)
    
    def derivative(self, x, h=10**(-5)):
        deriv = self.eval(x-2*h)-8*self.eval(x-h)+8*self.eval(x+h)-self.eval(x+2*h)
        return deriv/(12*h)