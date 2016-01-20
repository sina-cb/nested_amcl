#include "Observation.h"

size_t Observation::size(){
    return values.size();
}

Sample Observation::combine(Sample second){
    Sample result;

    for (size_t i = 0; i < this->size(); i++){
        result.values.push_back(values[i]);
    }
    for (size_t i = 0; i < second.size(); i++){
        result.values.push_back(second.values[i]);
    }

    return result;
}
