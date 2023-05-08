CasADi toolbox is used in this project. It is an open-source tool for nonlinear optimization and algorithmic differentiation. You can install it on https://web.casadi.org/get/.

Fmincon function in Matlab is also tried. But it is commented out for it took too long to run.
\bf{Main File}

"MPCPred.m".

\bf{GP Function files}

"GPRPredRQ.m"  GP Predict Function based on Rational quadratic kernel function. 

"GPRPredMatern.m"  GP Predict Function based on Matérn kernel function. You also need to assign aother parameter(1,3,5,7).

\bf{Fitrgp files}

"ObjFunc.m"

"NonLinearCons.m"
