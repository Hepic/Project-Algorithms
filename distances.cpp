#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "distances.h"

vector<double> find_closest_point(const vector<double> &curr_point, const vector<double> &lower_left_point, double delta) {
    int dim = curr_point.size();
    vector<double> closest_point;
    
    for (int i = 0; i < dim; ++i) {
        double val = curr_point[i] - lower_left_point[i] * delta;
        closest_point.push_back(lower_left_point[i]);
        
        if (val > delta / 2) {
            ++closest_point.back();
        }
    }
    
    return closest_point;
}

double euclidean_distance_square(const vector<double> &pnt_1, const vector<double> &pnt_2) {
    double dist = 0;

    for (int i = 0; i < (int)pnt_1.size(); ++i) {
        dist += (pnt_1[i] - pnt_2[i]) * (pnt_1[i] - pnt_2[i]);   
    }

    return dist;
}

vector<double> mean_point(const vector<double> &pnt_1, const vector<double> &pnt_2) {
    vector<double> mean;

    for (int i = 0; i < (int)pnt_1.size(); ++i) {
        double val = (pnt_1[i] + pnt_2[i]) / 2.0;
        mean.push_back(val);
    }

    return mean;
}

double discrete_frechet_distance(const Curve &curve_1, const Curve &curve_2, Curve &mean_traversal, bool path) {
    int N1 = curve_1.get_length();
    int N2 = curve_2.get_length();
    double **dp_solve = new double*[N1];
    double result, dist;
    
    for (int i = 0; i < N1; ++i) {
        dp_solve[i] = new double[N2];
    }
    
    dist = euclidean_distance_square(curve_1.get_point(0), curve_2.get_point(0));
    dp_solve[0][0] = dist;
    
    for (int i = 1; i < N1; ++i) {
        dist = euclidean_distance_square(curve_1.get_point(i), curve_2.get_point(0));
        dp_solve[i][0] = max(dp_solve[i - 1][0], dist);
    }
    
    for (int i = 1; i < N2; ++i) {
        dist = euclidean_distance_square(curve_1.get_point(0), curve_2.get_point(i));
        dp_solve[0][i] = max(dp_solve[0][i - 1], dist);
    }
    
    for (int i = 1; i < N1; ++i) {
        for (int j = 1; j < N2; ++j) {
            double dist = euclidean_distance_square(curve_1.get_point(i), curve_2.get_point(j));
            double val = -1;

            val = min(min(dp_solve[i - 1][j], dp_solve[i][j - 1]), dp_solve[i - 1][j - 1]);
            dp_solve[i][j] = max(val, dist);
        }
    }
    
    result = sqrt(dp_solve[N1 - 1][N2 - 1]);
    
    if (path) {
        int p1 = N1 - 1, p2 = N2 - 1;
        Curve reverse_mean_traversal;
        vector<double> mean;

        mean = mean_point(curve_1.get_point(p1), curve_2.get_point(p2));
        reverse_mean_traversal.insert_point(mean);
        
        while (p1 && p2) {
            if (dp_solve[p1 - 1][p2] < min(dp_solve[p1][p2 - 1], dp_solve[p1 - 1][p2 - 1])) {
                mean = mean_point(curve_1.get_point(--p1), curve_2.get_point(p2));
            } else if (dp_solve[p1][p2 - 1] < dp_solve[p1 - 1][p2 - 1]) {
                mean = mean_point(curve_1.get_point(p1), curve_2.get_point(--p2));
            } else {
                mean = mean_point(curve_1.get_point(--p1), curve_2.get_point(--p2));
            }
            
            reverse_mean_traversal.insert_point(mean); 
        }

        for (int i = (int)reverse_mean_traversal.get_length() - 1; i >= 0; --i) {
            mean_traversal.insert_point(reverse_mean_traversal.get_point(i));
        }
    }

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
    double dist = 0.0;
    
    if (!strcmp(dist_function, "DFT")) {
        dist = discrete_frechet_distance(curve_1, curve_2);
    } else if (!strcmp(dist_function, "DTW")) {
        dist = dynamic_time_wrapping(curve_1, curve_2);
    }
    
    return dist;
}
