import numpy as np
import matplotlib.pyplot as plt

# measurement: [value, -err, +err]
MAXITER = 1e5
verbose = False

#x = np.array([[7,2.323,2.989],
#     [3,1.416,2.080]])

#x = np.array([[9,2.676,3.342],
#     [1,0.6983,1.358]])
x = np.array([[9,2.676,3.342]])


def BarlowAsymmErrorAdder(x,modeltype="sigma",verbose=False,scaleErr=1.0):

    guessVal0th = 0.0
    guessErr0thPlus = 0.0
    guessErr0thMinus = 0.0
    if(modeltype not in ("sigma","variance")):
        print("Error! Modeltype should be one of 'sigma' or 'variance'")
        return

    for i in range(len(x)):
        guessVal0th += x[i][0]/len(x)
        guessErr0thPlus += 1/pow(x[i][2],2)
        guessErr0thMinus += 1/pow(x[i][1],2)

    guessErr0thPlus = np.sqrt(1/guessErr0thPlus)*scaleErr*np.random.rand(1)[0]
    guessErr0thMinus = np.sqrt(1/guessErr0thMinus)*scaleErr*np.random.rand(1)[0]
    print("Proceeding with initial guesses for mean: %g"%(guessVal0th))
    print("                                for +err: %g"%(guessErr0thPlus))
    print("                                for -err: %g"%(guessErr0thMinus))

    guessx = guessVal0th
    err = 1.0
    ctr = 0
    while err > 1e-6:
        L = 0.0
        numerator = 0
        denominator = 0
        for i in range(len(x)):
            if modeltype == "sigma":
                sigmai = 2*(x[i][2]*x[i][1])/(x[i][2]+x[i][1])
                sigmaprime_i = (x[i][2] - x[i][1])/(x[i][2] + x[i][1])
                wi = sigmai/pow(sigmai+sigmaprime_i*(guessx - x[i][0]),3.0)
                numerator += x[i][0]*wi
                denominator += wi
                L += -0.5*pow((guessx - x[i][0])/(sigmai+sigmaprime_i*(guessx - x[i][0])),2.0)
            elif modeltype == "variance":
                Vi = (x[i][2]*x[i][1])
                Vprime_i = (x[i][2] - x[i][1])
                wi = Vi/pow(Vi+Vprime_i*(guessx - x[i][0]),2.0)
                numerator += (x[i][0]-Vprime_i*pow(guessx-x[i][0],2.0)/(2.0*Vi))*wi
                denominator += wi
                L += -0.5*pow(guessx - x[i][0],2.0)/(Vi+Vprime_i*(guessx - x[i][0]))

        temp = guessx
        guessx = numerator/denominator
        err = abs(guessx-temp)
        if verbose:
            print("dx:",err,"lkhd:",L)
        ctr = ctr + 1
        if(ctr > MAXITER):
            print("Maximum iterations reached, breaking")
            print("Current value of dx:",err)
            break
    Lmax = L
    xmaxL = guessx
    print("MaxL:",xmaxL,Lmax,"Iter:",ctr)

    err = 1.0
    guessx = xmaxL + guessErr0thPlus
    dx = guessErr0thPlus
    ctr = 0
    while abs(err-0.5) > 1e-6:
        L = 0.0
        for i in range(len(x)):
            if modeltype == "sigma":
                sigmai = 2*(x[i][2]*x[i][1])/(x[i][2]+x[i][1])
                sigmaprime_i = (x[i][2] - x[i][1])/(x[i][2] + x[i][1])
                L += -0.5*pow((guessx - x[i][0])/(sigmai+sigmaprime_i*(guessx - x[i][0])),2.0)
            elif modeltype == "variance":
                Vi = (x[i][2]*x[i][1])
                Vprime_i = (x[i][2] - x[i][1])
                L += -0.5*pow(guessx - x[i][0],2.0)/(Vi+Vprime_i*(guessx - x[i][0]))

        preverr = err
        err = (Lmax-L)
        if verbose:
            print("Err:",err,"Lmax:",Lmax,"L:",L,"x:",guessx,"dx:",dx)

        if(err > 0.5):
            guessx -= dx
            if(preverr < 0.5):
                dx = dx/10
        else:
            guessx = guessx + dx        
        ctr = ctr + 1
        if(ctr > MAXITER):
            print("Maximum iterations reached, breaking")
            print("Current value of dL:",err)
            break

    xUp = guessx
    print("UpError:",xUp - xmaxL,"Iter:",ctr,"dL:",err)

    err = 1.0
    guessx = xmaxL - guessErr0thMinus
    dx = guessErr0thMinus
    ctr = 0
    while abs(err-0.5) > 1e-6:
        L = 0.0
        for i in range(len(x)):
            if modeltype == "sigma":
                sigmai = 2*(x[i][2]*x[i][1])/(x[i][2]+x[i][1])
                sigmaprime_i = (x[i][2] - x[i][1])/(x[i][2] + x[i][1])
                L += -0.5*pow((guessx - x[i][0])/(sigmai+sigmaprime_i*(guessx - x[i][0])),2.0)
            elif modeltype == "variance":
                Vi = (x[i][2]*x[i][1])
                Vprime_i = (x[i][2] - x[i][1])
                L += -0.5*pow(guessx - x[i][0],2.0)/(Vi+Vprime_i*(guessx - x[i][0]))

        preverr = err
        err = (Lmax-L)
        if verbose:
             print("Err:",err,"Lmax:",Lmax,"L:",L,"guessx:",guessx,"dx:",dx)

        if(err > 0.5): #while we go to the negative -0.5, we need a zero crossing for us to change the iterator size
            guessx += dx
            if(preverr < 0.5):
                dx = dx/10
        else:
            guessx = guessx - dx        
        ctr = ctr + 1
        if(ctr > MAXITER):
            print("Maximum iterations reached, breaking")
            print("Current value of dL:",err)
            break

    xDown = guessx
    print("DownError:",xDown - xmaxL,"Iters:",ctr,"dL:",err)
    return([xmaxL,xDown-xmaxL,xUp-xmaxL],Lmax)

def main():
    one = BarlowAsymmErrorAdder(x,"sigma",False,scaleErr=0.033243)
    two = BarlowAsymmErrorAdder(x,"variance",True,scaleErr=0.033243)
    print(one, two)

if __name__ == "__main__":
    main()
