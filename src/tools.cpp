#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    // ... your code here
    if (estimations.size() == 0)
    {
        cout << "estimation vector size is zero";
    }else if(estimations.size() != ground_truth.size())
    {
        cout << "estimation vector size is not equal to ground truth vector size";
    }

    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        VectorXd est = estimations[i];
        VectorXd gt = ground_truth[i];

        VectorXd a = est - gt;
        VectorXd b = a.array() * a.array();

        rmse += b;
    }

    //calculate the mean
    rmse = rmse / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}
