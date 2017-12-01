#ifndef CURVE_H
#define CURVE_H

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class Curve {
    int dim, int_id;
    string id;
    vector<vector<double> > curve;

    public:
        Curve();
        Curve(string, int, int);
        Curve(string, int, int, const vector<vector<double> >&);
        void set_id(string, int);
        void insert_point(const vector<double>&);
        void clear_curve();
        int get_dimension() const;
        int get_length() const;
        string get_id() const;
        int get_int_id() const;
        double get_coord_point(int, int) const;
        const vector<double>& get_point(int) const;
        const vector<double>& get_last_point() const;
        vector<double> get_convert_vector() const;
        bool is_empty() const;
        void print_curve() const;
        void append_curve(const Curve&);
        bool equal_curves(const Curve&) const;
        bool operator<(const Curve &curve) const;
};

extern vector<Curve> input_curves;

#endif
