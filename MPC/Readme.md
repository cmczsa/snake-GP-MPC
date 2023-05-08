# CasADi toolbox

CasADi toolbox is used in this project. It is an open-source tool for nonlinear optimization and algorithmic differentiation. You can install it on https://web.casadi.org/get/.

Fmincon function in Matlab is also tried. But it is commented out for it took too long to run.

# Files

__1. Main File__

"MPCPred.m".

__2. GP Function Files__

"GPRPredRQ.m"  GP Predict Function based on Rational quadratic kernel function. 

"GPRPredMatern.m"  GP Predict Function based on Matérn kernel function. You also need to assign aother parameter(1,3,5,7).

__3. Model Files__

"MPC_CpxModelling.m"

Complicated snake robot model file which constructed in folder "ComplicatedModel".


__4. Fitrgp Files__

"ObjFunc.m"

"NonLinearCons.m"

__5. Other Files__

.mat and .xlsx files are the trained models and the computed data.

# More MPC Examples

If you want to learn more examples about casADi and MPC, please browse my another repository "".
