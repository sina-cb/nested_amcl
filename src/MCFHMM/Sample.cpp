#include "Sample.h"
#include <random>
#include <chrono>

Sample::~Sample(){
    values.clear();
//    LOG(INFO) << "SAMPLE CLASS DESTRUCTOR!";
}

size_t Sample::size(){
    return values.size();
}

void Sample::init_rand(vector<double> *low_limit, vector<double> *high_limit){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    for (size_t i = 0; i < low_limit->size(); i++){
        uniform_real_distribution<double> dist((*low_limit)[i], (*high_limit)[i]);
        values.push_back(dist(gen));
    }
}

Sample Sample::combine(vector<double> second){
    Sample result;

    for (size_t i = 0; i < this->size(); i++){
        result.values.push_back(values[i]);
    }
    for (size_t i = 0; i < second.size(); i++){
        result.values.push_back(second[i]);
    }

    return result;
}

string Sample::str(){
    stringbuf buf;
    ostream os (&buf);
    os << "\tValues -> ";
    for (size_t i = 0; i < values.size(); i++){
        os << values[i] << " ";
    }
    os << "\tP -> " << p << endl;

    return buf.str();
}
