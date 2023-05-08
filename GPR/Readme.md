# Introduction

This folder aims to train and predict the gait pattern of a snake robot. 

# GPML Toolbox

GPML toolbox is used in this project. You can install it on  http://gaussianprocess.org/gpml/code/matlab/doc/.

We also used fitrgp function in Matlab, while it can't be customed. 

# Files

__1. Main File__

"GPR.m".

__2. GP Function Files__

"GPRPredRQ.m" 

GP Predict Function based on Rational quadratic kernel function.

"GPRPredMatern.m" 

GP Predict Function based on Matérn kernel function. You also need to assign aother         parameter(1,3,5,7).

__3. Model Files__

"MPC_CpxModelling.m" 

Complicated snake robot model file which constructed in folder "ComplicatedModel".

__4. Other Files__

.mat and .xlsx files are the trained models and the computed data.
