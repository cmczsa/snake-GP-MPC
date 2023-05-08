This folder aims to train and predict the gait pattern of a snake robot. 

GPML toolbox is used in this project. You can install it on  http://gaussianprocess.org/gpml/code/matlab/doc/.

We also used fitrgp function in Matlab, while it can't be customed. 

__Main File__

"GPR.m".

__GP Function Files__

"GPRPredRQ.m" 

GP Predict Function based on Rational quadratic kernel function.

"GPRPredMatern.m" 

GP Predict Function based on Matérn kernel function. You also need to assign aother         parameter(1,3,5,7).

__Nodel Files__

"MPC_CpxModelling.m" 

Complicated snake robot model file which constructed in folder "ComplicatedModel".

__Other Files__

.mat and .xlsx files are the trained models and the computed data.
