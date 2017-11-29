#include "hash_functions.h"
#include "distances.h"
#include "help_functions.h"

vector<int> vec_r;
vector<vector<vector<double> > > t_shift, vec_line;

vector<Curve> grid_hashing(double delta, const vector<double> &t_shift_grid) {
    vector<Curve> curves = input_curves;
    vector<Curve> modified_curves;

    for (int i = 0; i < curves.size(); ++i) {
        Curve image_curve;

        for (int j = 0; j < curves[i].get_length(); ++j) {
            vector<double> lower_left_point, shifted_point;
            
            for(int k = 0; k < curves[i].get_dimension(); ++k) {
                double value = curves[i].get_coord_point(k, j);
                double shifted_value = value - t_shift_grid[k];
                double image_value = (int)(shifted_value / delta);
                
                lower_left_point.push_back(image_value);
                shifted_point.push_back(shifted_value);
            }
            
            vector<double> closest_point = find_closest_point(shifted_point, lower_left_point, delta);

            if (!image_curve.is_empty() && image_curve.get_last_point() == closest_point) { // remove duplicates
                continue;
            }
            
            image_curve.insert_point(closest_point);
        }
        
        modified_curves.push_back(image_curve);
    }
    
    return modified_curves;
}

Curve grid_hashing_curve(double delta, const vector<double> &t_shift_grid, const Curve &curve) {
    Curve modified_curve;
    modified_curve.set_id(curve.get_id());
    
    for (int i = 0; i < curve.get_length(); ++i) {
        vector<double> lower_left_point, shifted_point;
        
        for(int j = 0; j < curve.get_dimension(); ++j) {
            double value = curve.get_coord_point(j, i);
            double shifted_value = value - t_shift_grid[j];
            double image_value = (int)(shifted_value / delta);
            
            lower_left_point.push_back(image_value);
            shifted_point.push_back(shifted_value);
        }
        
        vector<double> closest_point = find_closest_point(shifted_point, lower_left_point, delta);

        if (!modified_curve.is_empty() && modified_curve.get_last_point() == closest_point) { // remove duplicates
            continue;
        }
        
        modified_curve.insert_point(closest_point);
    }
    
    return modified_curve;
}

vector<Curve> multiple_grids(double delta, int k, int hash_id) {
    vector<Curve> curves = input_curves;
    vector<Curve> concat_curves(curves.size(), Curve()), modified_curves;
    vector<double> t_shift_grid;
    vector<vector<double> > t_shift_k_grids;
    int dim;
   
    if (!curves.empty()) {
        dim = curves[0].get_dimension();
    } 

    for (int i = 0; i < k; ++i) {
        t_shift_grid.clear();
    
        for (int j = 0; j < dim; ++j) {
            double val = uniform_distribution(0, dim);
            t_shift_grid.push_back(val);
        }

        t_shift_k_grids.push_back(t_shift_grid);
        modified_curves = grid_hashing(delta, t_shift_grid);

        for (int j = 0; j < modified_curves.size(); ++j) {
            concat_curves[j].append_curve(modified_curves[j]);
        }
    } 

    t_shift.push_back(t_shift_k_grids);

    for (int i = 0; i < curves.size(); ++i) {
        concat_curves[i].set_id(curves[i].get_id());
    }

    return concat_curves;
}

Curve multiple_grids_curve(double delta, int k, int hash_id, const Curve &curve) {
    Curve concat_curve;
    concat_curve.set_id(curve.get_id());

    for (int i = 0; i < k; ++i) { // create concated grid_curve after k grid_hashings
        Curve modified_curve = grid_hashing_curve(delta, t_shift[hash_id][i], curve);
        concat_curve.append_curve(modified_curve);
    }

    return concat_curve;
}

vector<double> prob_lsh_euclid_hashing(const vector<double> &vec, int w, int k_vec, int hash_id) {
    vector<double> vec_lsh;

    for (int i = 0; i < k_vec; ++i) {
        double dot = dot_product(vec, vec_line[hash_id][i]);
        double t = uniform_distribution(0, w);

        int val = (dot + t) / w;
        vec_lsh.push_back(val);
    }
    
    return vec_lsh;
}

long long vector_hashing(const vector<double> &vec) {
    long long sum = 0, MOD = 4294967291LL;

    for (int i = 0; i < vec.size(); ++i) {
        if (i >= vec_r.size()) {
            vec_r.push_back(rand() % MAX_R);
        }

        sum += (((long long)vec_r[i] + MOD) % MOD) * (((long long)vec[i] + MOD) % MOD);
        sum %= MOD;
    }
    
    sum = (sum + MOD) % MOD;
    return sum;
}
