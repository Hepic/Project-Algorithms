#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "distances.h"

double euclidean_distance_square(const vector<double> &pnt_1, const vector<double> &pnt_2) {
    double dist = 0;

    for (int i = 0; i < pnt_1.size(); ++i) {
        dist += (pnt_1[i] - pnt_2[i]) * (pnt_1[i] - pnt_2[i]);   
    }

    return dist;
}

double discrete_frechet_distance(const Curve &curve_1, const Curve &curve_2) {
    int N1 = curve_1.get_length();
    int N2 = curve_2.get_length();
    double **dp_solve = new double*[N1];
    double result;
    
    for (int i = 0; i < N1; ++i) {
        dp_solve[i] = new double[N2];
        
        for (int j = 0; j < N2; ++j) {
            dp_solve[i][j] = 0;
        }
    }
     
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            double dist = euclidean_distance_square(curve_1.get_point(i), curve_2.get_point(j));
            double val = -1;

            if (i-1 >= 0) {
                val = dp_solve[i-1][j];
            }
            
            if (j-1 >= 0) {
                val = (val == -1 ? dp_solve[i][j-1] : min(val, dp_solve[i][j-1]));
            }
            
            if (i-1 >= 0 && j-1 >= 0) {
                val = (val == -1 ? dp_solve[i-1][j-1] : min(val, dp_solve[i-1][j-1]));
            }
            
            dp_solve[i][j] = max(val, dist);
        }
    }
    
    result = sqrt(dp_solve[N1-1][N2-1]);
    
    for (int i = 0; i < N1; ++i) {
        delete[] dp_solve[i];
    }

    delete[] dp_solve;
    return result;
}

double dynamic_time_wrapping(const Curve &curve_1, const Curve &curve_2) {
    int N1 = curve_1.get_length();
    int N2 = curve_2.get_length();
    double **dp_solve = new double*[N1];
    double result;

    for (int i = 0; i < N1; ++i) {
        dp_solve[i] = new double[N2];

        for (int j = 0; j < N2; ++j) {
            dp_solve[i][j] = 0;
        }
    }
    
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            double dist = euclidean_distance_square(curve_1.get_point(i), curve_2.get_point(j));
            double val = -1;

            if (i-1 >= 0) {
                val = dp_solve[i-1][j];
            }

            if (j-1 >= 0) {
                val = (val == -1 ? dp_solve[i][j-1] : min(val, dp_solve[i][j-1]));
            }
            
            if (i-1 >= 0 && j-1 >= 0) {
                val = (val == -1 ? dp_solve[i-1][j-1] : min(val, dp_solve[i-1][j-1]));
            }
           
            dp_solve[i][j] = dist;
            
            if (val != -1) {
                dp_solve[i][j] += val;
            }
        }
    }
    
    result = sqrt(dp_solve[N1-1][N2-1]);
    
    for (int i = 0; i < N1; ++i) {
        delete[] dp_solve[i];
    }

    delete[] dp_solve;    
    return result;
}

double compute_distance(const Curve &curve_1, const Curve &curve_2, const char *dist_function) {
    double dist;
    
    if (!strcmp(dist_function, "DFT")) {
        dist = discrete_frechet_distance(curve_1, curve_2);
    } else if (!strcmp(dist_function, "DTW")) {
        dist = dynamic_time_wrapping(curve_1, curve_2);
    }
    
    return dist;
}
