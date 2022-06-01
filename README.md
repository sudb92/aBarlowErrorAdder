aBarlowErrorAdder
----------------

* Very quick NumPy implementation of R Barlow's method(https://arxiv.org/abs/physics/0406120) to combine many measurements with asymmetric statistical errorbars
* Implements both the 'linear standard deviation' and 'variance' methods described in the original reference
* Put the file in your project directory, and use the line 'import BarlowErrorModule' or equivalents to continue
* BarlowErrorModule.BarlowAsymmErrorAdder() is the function call - it takes the following arguments
	- "x" : a NumPy array that has the form [[meas0,err0-,err0+], [meas1,err1-,err1+], [meas2,err2-,err2+], ... ,[measN,errN-,errN+]]
	 		(Here, meas0 is the central value, err0- is the -ve errorbar, err0+ is its +ve errorbar) 
	- "modeltype": Can be  "sigma"(default) or "variance". Chooses the two models described in the original paper.
	- "verbose" : True/False which decides if the iteration progress is printed out or otherwise
	- "scaleErr": It's useful in some cases to override the initial guess size when calculating the dlog(L) = -0.5 points. Since there is a built-in
				  randomization to the stepped value, this can be used to fix its order of magnitude. Useful to give it a smaller number if the
				  first guess overshoots the dlog(L) = -0.5 point, a larger value if the iterations take too long. It's a good idea to look at the verbose output
				  before proceeding. 
				  
All implemented capabilities are almost exactly as described in the original reference. The script can also be run stand-alone - as 

python3 BarlowErrorModule.py

which gives results from a test case.

Any feedback, questions, comments, criticisms are welcome - and appreciated! 


---
B. Sudarsan <br>
sbalak2@lsu.edu <br>
bsudarsan92@gmail.com <br>
LSU Baton Rouge
