#-----------------------------------------------------------------------------------------------
class Function:
    def __init__(self):
        pass
    def df(self, x):
        return 2.0*x+10**(-20)

    def f(self, x, f0):
        f = x**2-f0
        return f

#-----------------------------------------------------------------------------------------------

class RootFinder:
    def __init__(self):
        pass

    #-----------------------------------------------------------------------------------------------

    def NewtonRalphson(self,f, df,f0, x0):
        x_last = x0
        run = True

        i=0

        while run:
            i += 1
            xi = x_last - f(x_last, f0)/df(x_last)

            if abs ((xi-x_last)/(x_last+10e-15)) <= 10**-15:
                run = False
                return [i, xi]

            x_last = xi

    #-----------------------------------------------------------------------------------------------

    def Bisection(self, f,f0, a, b):

        if abs(f(a, f0)) <= 10e-15:
            return a

        if abs(f(b, f0)) <= 10e-5:
            return b

        if f(a, f0)*f(b, f0) > 0:
            print("Solution for f = {:5.4e} is not within [{}, {}]. Try new initial values. ".format(f0,a,b))
            exit(1)


        else:
            i = 0
            run = True
            c_last = 0.0

            while run:

                i += 1

                c = (a+b)/2.0

                if abs((c-c_last)/c) <=10e-20:
                    return c

                elif f(c, f0)*f(a, f0)<0:
                    b = c

                elif f(c, f0)*f(b, f0)<0:
                    a = c

                c_last = c

            return c

#-----------------------------------------------------------------------------------------------


def main():
    rootFinderObj = RootFinder()
    functionObj = Function()
    f0 = 100.0
    #root = rootFinderObj.NewtonRalphson(functionObj.f, functionObj.df, -0.1)
    root = rootFinderObj.Bisection(functionObj.f,f0, 0.0, 11.0)

    print ("iterations: {}\t root: {}".format(root[0], root[1]))

#-----------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
