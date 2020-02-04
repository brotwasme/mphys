#check the workings off the class

class over:
    def __init__(self, pop="pop"):
        self.pop = pop
        self.a = 2
        self.b = 1

    def over_method(self,d=1):
        return self.a + self.b, d

class under(over):
    def __init__(self, pat="pat"):
        super(under,self).__init__(pop="yay")
        self.c = 4
    
    def over_method(self):
        outs, d = super(under,self).over_method()
        return outs, 2, d, "done"

a = under()
print(a.over_method())