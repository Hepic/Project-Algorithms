#include "curve.h"

vector<Curve> input_curves;

Curve::Curve() {
    this->id = "";
    this->int_id = -1;
}

Curve::Curve(string id, int int_id, int dim) {
    this->id = id;
    this->int_id = int_id;
    this->dim = dim;
}

Curve::Curve(string id, int int_id, int dim,const vector<vector<double> > &curve) {
    this->id = id;
    this->int_id = int_id;
    this->dim = dim;
    this->curve = curve;
}

void Curve::set_id(string id, int int_id) {
    this->id = id;
    this->int_id = int_id;
}

void Curve::set_dimension(int dim) {
    this->dim = dim;
}

void Curve::insert_point(const vector<double> &point) {
    curve.push_back(point);
}

void Curve::clear_curve() {
    curve.clear();
}

int Curve::get_dimension() const {
    return dim;
}

int Curve::get_length() const {
    return curve.size();
}

string Curve::get_id() const {
    return id;
}

int Curve::get_int_id() const {
    return int_id;
}

double Curve::get_coord_point(int coord, int pnt) const {
    return curve[pnt][coord];
}

const vector<double>& Curve::get_point(int pnt) const {
    return curve[pnt];
}

const vector<double>& Curve::get_last_point() const {
    return curve.back();
}

vector<double> Curve::get_convert_vector() const {
    vector<double> gen_vector;
    
    for (int i = 0; i < (int)curve.size(); ++i) {
        for (int j = 0; j < (int)curve[i].size(); ++j) {
            gen_vector.push_back(curve[i][j]);
        }
    }
    
    return gen_vector;
}

bool Curve::is_empty() const {
    return curve.empty();
}

void Curve::print_curve() const {
    cout << "Id: " << id << "\n";
    cout << "int_Id: " << int_id << "\n";

    for (int i = 0; i < (int)curve.size(); ++i) {
        cout << "Point " << i << ": ";

        for (int j = 0; j < (int)curve[i].size(); ++j) {
            cout << curve[i][j] << " ";
        }

        cout << "\n";
    }

    cout << "\n";
}

void Curve::append_curve(const Curve &new_curve) {
    for (int i = 0; i < new_curve.get_length(); ++i) {
        curve.push_back(new_curve.get_point(i));
    }
}

bool Curve::equal_curves(const Curve &curve_1) const {
    if ((int)curve.size() == curve_1.get_length()) {
        for (int i = 0; i < (int)curve.size(); ++i) {
            if (curve[i] != curve_1.get_point(i)) {
                return false;
            }
        }
        
        return true;
    }
    
    return false;
}

bool Curve::operator<(const Curve &curve) const {
    return this->id < curve.get_id();
}
