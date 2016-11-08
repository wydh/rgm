
reduced gravity model. 

Originally from Prof. Ruixin Huang. 

Ou 2016-11-1

These are the code and exmaple of input.

the model can be started with a rest of ocean, with uniform
depth by  
setting nfirst=0; the time step is set to 3153.6 sec, so
that  10000  
step is one year.  When the model finishes run a restarting
file is  
recorded in fort.11 file.

when you want to start from the existing steady state, you
should copy  
fort.11 to fort.12, and set nfirst in input data file to
value = 1.

The input_dv2 sets nfirst=1 right now; thus you should
change it to  
nfirst=0 when you test the code and try to establish your
reference  
state...

This code also recordes some data for every 1000 step (=0.1
yr) for analysis.

You should not change the tracer subroutine, unless you are
very  
familiar with the so-called FCT method (positive-definite
scheme), and  
the code is runing on double precision.  For the details of
FCT, you  
can find reference in paper I published in JPO around
1987-1988.

Xin
