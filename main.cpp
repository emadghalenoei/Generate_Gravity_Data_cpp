//
//  main.cpp
//
//  Created by Emad Ghalenoei in 2020.
//  Copyright (c) 2020 Emad.Ghalenoei. All rights reserved.



#include <fstream> // for file access
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <True_Modelath> // For Absolute function
//#include "matrix.h" // Matrix Class
#include <random>
#include <chrono>
#include "matrix.hpp" // Matrix Class

#define PI 3.14159265

using namespace std;

int main()
{


    Matrix profile=Matrix("2D_profile.txt");  // data profile
    Matrix kernel=Matrix("kernel.txt"); // Kernel (You need to define kernel or load it)
    Matrix True_Model=Matrix("True_Model.txt"); //True Model

    int num_obs = profile.getRows();
    Matrix x_s = profile.MatExtraction(1, profile.getRows(), 1, 1);
    Matrix x=x_s;
    double Z_model_space = 100;      // difference betweeen each depth layer (m)
    double Z1=0;                     // first depth at model space (m)
	double Z2=10000;                 // final depth at model space (m)
    //int num_layer=100;                 // number of z layer
    //Matrix z (num_layer,1);
    //z(0,0)=Z1+(Z_model_space/2);
    //z = Matrix::linespace(z(0,0), Z_model_space, num_layer);
	Matrix z = Matrix::colon(Z1+(Z_model_space/2),Z_model_space,Z2);
    Matrix X(z.getRows(),x.getRows());
    Matrix Z(z.getRows(),x.getRows());
    tie(X, Z) = Matrix::meshgrid(x, z);
    double density = -0.4;   // Unit(kg/m^3)
    density = density*1000;  // Unit(kg/m^3)
    True_Model=True_Model * density;
    //int tot_parameters= True_Model.getLength();
    Matrix True_Model_vec=True_Model.vectorize ();
    Matrix dg = kernel * True_Model_vec;
    dg = dg * 100000;      // 1 SI = 1e5 mGal



    Matrix ind = Matrix::colon(0,1,num_obs-1);
    Matrix Dij=ind.toeplitz();
    double re=15;                  // lower this number leads to converging to 0 at right side of function
    double L0=20;                          // it brings the function into lower number (vertically), higher L1, upper number
    Matrix T1 = (Dij*(2*PI/L0)).cos_mat();
    Matrix T2 =(Dij/(-1*re)).exp_mat();
    Matrix T = Matrix::dot(T1,T2);
    double noise_level = 0.03;
    Matrix sigma = dg.absolute() * noise_level;
    Matrix T3 =sigma.transpose();
    Matrix T4 = sigma * T3;
    Matrix C = Matrix::dot(T,T4);
    C = C.damping();
    Matrix mean(num_obs,1);
    Matrix std_mat = mean +1;
    Matrix np = Matrix::normrnd(mean, std_mat);
    Matrix L = C.cholesky("Lower");
    Matrix n = L * np;
    Matrix dg_obs=dg+n;
    Matrix noise_lvl (num_obs,1);
    noise_lvl = noise_lvl.one();
    noise_lvl = noise_lvl * noise_level;
    Matrix table = Matrix::join_in_right({x_s, dg_obs, dg,noise_lvl });
    remove ("table.txt");
    table.write2txtfile("table.txt");

    remove ("C_original.txt");
    C.write2txtfile("C_original.txt");

    return 0;

}
