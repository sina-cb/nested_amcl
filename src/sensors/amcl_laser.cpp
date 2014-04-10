/*
 *  Player - One Hell of a Robot Server
 *  Copyright (C) 2000  Brian Gerkey   &  Kasper Stoy
 *                      gerkey@usc.edu    kaspers@robotics.usc.edu
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
///////////////////////////////////////////////////////////////////////////
//
// Desc: AMCL laser routines
// Author: Andrew Howard
// Date: 6 Feb 2003
// CVS: $Id: amcl_laser.cc 7057 2008-10-02 00:44:06Z gbiggs $
//
///////////////////////////////////////////////////////////////////////////

#include <sys/types.h> // required by Darwin
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include "amcl_laser.h"

#include <ros/console.h>

#define LASER_WEIGHTAGE 1
#define COLOR_WEIGHTAGE 0

using namespace amcl;

////////////////////////////////////////////////////////////////////////////////
// Default constructor
AMCLLaser::AMCLLaser(size_t max_beams, map_t* map) : AMCLSensor()
{
    this->time = 0.0;

    this->max_beams = max_beams;
    this->map = map;

    return;
}


// Default constructor with color_map added
AMCLLaser::AMCLLaser(size_t max_beams, map_t* map, map_t* color_map) : AMCLSensor()
{
    this->time = 0.0;

    this->max_beams = max_beams;
    this->map = map;
    this->color_map = color_map;

    return;
}




void
AMCLLaser::SetModelBeam(double z_hit,
                        double z_short,
                        double z_max,
                        double z_rand,
                        double sigma_hit,
                        double lambda_short,
                        double chi_outlier)
{
    this->model_type = LASER_MODEL_BEAM;
    this->z_hit = z_hit;
    this->z_short = z_short;
    this->z_max = z_max;
    this->z_rand = z_rand;
    this->sigma_hit = sigma_hit;
    this->lambda_short = lambda_short;
    this->chi_outlier = chi_outlier;
}

void
AMCLLaser::SetModelLikelihoodField(double z_hit,
                                   double z_rand,
                                   double sigma_hit,
                                   double max_occ_dist)
{
    this->model_type = LASER_MODEL_LIKELIHOOD_FIELD;
    this->z_hit = z_hit;
    this->z_max = z_max;
    this->z_rand = z_rand;
    this->sigma_hit = sigma_hit;

    map_update_cspace(this->map, max_occ_dist);
    map_update_cspace(this->color_map, max_occ_dist);
}


////////////////////////////////////////////////////////////////////////////////
// Apply the laser sensor model
bool AMCLLaser::UpdateSensor(pf_t *pf, AMCLSensorData *data)
{
    AMCLLaserData *ndata;

    ndata = (AMCLLaserData*) data;
    if (this->max_beams < 2)
        return false;

    // Apply the laser sensor model
    if(pf->nesting_lvl < 1){
        if(this->model_type == LASER_MODEL_BEAM)
            pf_update_sensor(pf, (pf_sensor_model_fn_t) BeamModel, data);
        else if(this->model_type == LASER_MODEL_LIKELIHOOD_FIELD)
            pf_update_sensor(pf, (pf_sensor_model_fn_t) LikelihoodFieldModel, data);
        else
            pf_update_sensor(pf, (pf_sensor_model_fn_t) BeamModel, data);
    }

    else{
        //ROS_INFO("inside nested UpdateSensor");
        if(this->model_type == LASER_MODEL_BEAM)
            pf_update_nested_sensor(pf, (pf_sensor_model_fn_t) BeamModel, (pf_nested_sensor_model_fn_t) NestedBeamModel, data);
        else if(this->model_type == LASER_MODEL_LIKELIHOOD_FIELD)
            pf_update_nested_sensor(pf, (pf_sensor_model_fn_t) LikelihoodFieldModel, (pf_nested_sensor_model_fn_t) NestedBeamModel, data);
        else
            pf_update_nested_sensor(pf, (pf_sensor_model_fn_t) BeamModel, (pf_nested_sensor_model_fn_t) NestedBeamModel, data);
    }
    /*
    if(pf->nesting_lvl > 0){
        for(int i=0; i < pf->sets[pf->current_set].sample_count ; i++){
            pf_t *nested_pf_set, *nested_pf;
            nested_pf_set = pf_get_this_nested_set(pf, pf->current_set);
            nested_pf = nested_pf_set + i;
            pf_update_sensor(nested_pf, (pf_sensor_model_fn_t) NestedBeamModel, data);
        }
    }
*/
    return true;
}


////////////////////////////////////////////////////////////////////////////////
// Determine the probability for the given pose
double AMCLLaser::BeamModel(AMCLLaserData *data, pf_sample_set_t* set)
{
    AMCLLaser *self;
    int i, j, step;
    double z, pz;
    double p;
    double map_range;
    double obs_range, obs_bearing;
    double total_weight;
    pf_sample_t *sample;
    pf_vector_t pose;

    self = (AMCLLaser*) data->sensor;

    total_weight = 0.0;

    // Compute the sample weights
    for (j = 0; j < set->sample_count; j++)
    {
        sample = set->samples + j;
        pose = sample->pose;

        // Take account of the laser pose relative to the robot
        pose = pf_vector_coord_add(self->laser_pose, pose);

        p = 1.0;

        step = (data->range_count - 1) / (self->max_beams - 1);
        for (i = 0; i < data->range_count; i += step)
        {
            obs_range = data->ranges[i][0];
            obs_bearing = data->ranges[i][1];


            //@KPM TODO : add color testing to map_calc_range here

            // Compute the range according to the map
            map_range = map_calc_range(self->map, pose.v[0], pose.v[1],
                                       pose.v[2] + obs_bearing, data->range_max);
            pz = 0.0;

            // Part 1: good, but noisy, hit
            z = obs_range - map_range;
            pz += self->z_hit * exp(-(z * z) / (2 * self->sigma_hit * self->sigma_hit));

            // Part 2: short reading from unexpected obstacle (e.g., a person)
            if(z < 0)
                pz += self->z_short * self->lambda_short * exp(-self->lambda_short*obs_range);

            // Part 3: Failure to detect obstacle, reported as max-range
            if(obs_range == data->range_max)
                pz += self->z_max * 1.0;

            // Part 4: Random measurements
            if(obs_range < data->range_max)
                pz += self->z_rand * 1.0/data->range_max;

            // TODO: outlier rejection for short readings

            assert(pz <= 1.0);
            assert(pz >= 0.0);
            //      p *= pz;
            // here we have an ad-hoc weighting scheme for combining beam probs
            // works well, though...
            p += pz*pz*pz;
        }

        sample->weight *= p;
        total_weight += sample->weight;
    }

    return(total_weight);
}

double AMCLLaser::LikelihoodFieldModel(AMCLLaserData *data, pf_sample_set_t* set)
{
    AMCLLaser *self;
    int i, j, step;
    double z, pz;
    double color_z; //, color_pz; //KPM: using the same pz for both (combining earlier itself)
    double p;
    double obs_range, obs_bearing, obs_color;
    double total_weight;
    pf_sample_t *sample;
    pf_vector_t pose;
    pf_vector_t hit;

    self = (AMCLLaser*) data->sensor;

    total_weight = 0.0;

    // Compute the sample weights
    for (j = 0; j < set->sample_count; j++)
    {
        sample = set->samples + j;
        pose = sample->pose;

        // Take account of the laser pose relative to the robot
        pose = pf_vector_coord_add(self->laser_pose, pose);

        p = 1.0;

        // Pre-compute a couple of things
        double z_hit_denom = 2 * self->sigma_hit * self->sigma_hit;
        double z_rand_mult = 1.0/data->range_max;

        step = (data->range_count - 1) / (self->max_beams - 1);
        for (i = 0; i < data->range_count; i += step)
        {
            obs_range = data->ranges[i][0];
            obs_bearing = data->ranges[i][1];
            obs_color = data->ranges[i][2];

            // This model ignores max range readings
            if(obs_range >= data->range_max)
                continue;

            pz = 0.0;

            // Compute the endpoint of the beam
            hit.v[0] = pose.v[0] + obs_range * cos(pose.v[2] + obs_bearing);
            hit.v[1] = pose.v[1] + obs_range * sin(pose.v[2] + obs_bearing);

            // Convert to map grid coords.
            int mi, mj;
            mi = MAP_GXWX(self->map, hit.v[0]);
            mj = MAP_GYWY(self->map, hit.v[1]);

            // Part 1: Get distance from the hit to closest obstacle.
            // Off-map penalized as max distance
            if( ( !MAP_VALID(self->map, mi, mj) || (self->map->cells[MAP_INDEX(self->map,mi,mj)].occ_state > -1) )){
                z = self->map->max_occ_dist;
                color_z = self->color_map->max_occ_dist;
            }
            else{
                z = self->map->cells[MAP_INDEX(self->map,mi,mj)].occ_dist;
                color_z = self->color_map->cells[MAP_INDEX(self->color_map, mi,mj)].occ_dist;
            }

            if(obs_color == 50.0){   // KPM: adding checking for color matches

                // Gaussian model
                // NOTE: this should have a normalization of 1/(sqrt(2pi)*sigma)
                if(LASER_WEIGHTAGE == 1){
                    pz = pz + (self->z_hit * exp(-(z * z) / z_hit_denom));
                }
                if(COLOR_WEIGHTAGE == 1){
                    pz *= (self->z_hit * exp(-(color_z * color_z) / z_hit_denom));
                }


                // Part 2: random measurements
                pz = pz + (self->z_rand * z_rand_mult) *(self->z_rand * z_rand_mult);

            }

            else{   // KPM: adding checking for color matches

                // Gaussian model
                // NOTE: this should have a normalization of 1/(sqrt(2pi)*sigma)
                if(LASER_WEIGHTAGE == 1){
                    pz = pz + (self->z_hit * exp(-(z * z) / z_hit_denom));
                }
                if(COLOR_WEIGHTAGE == 1){
                    pz *= (1- (self->z_hit * exp(-(color_z * color_z) / z_hit_denom)));
                }


                // Part 2: random measurements
                pz = pz + (self->z_rand * z_rand_mult) *(self->z_rand * z_rand_mult);
            }


            // TODO: outlier rejection for short readings

            assert(pz <= 1.0);
            assert(pz >= 0.0);
            //      p *= pz;
            // here we have an ad-hoc weighting scheme for combining beam probs
            // works well, though...
            p += pz*pz*pz;
        }

        sample->weight *= p;
        //ROS_INFO("sample weight: %f", sample->weight);
        total_weight += sample->weight;
    }

    return(total_weight);
}




////////////////////////////////////////////////////////////////////////////////
// Determine the probability of nested particles for the given pose
double AMCLLaser::NestedBeamModel(pf_sample_t *upper_sample, AMCLLaserData *data, pf_sample_set_t* set)
{
    AMCLLaser *self;
    int i, j, step;
    double z, pz, pz_color;
    double p;
    //    double map_range;
    double obs_range, obs_bearing, obs_color;
    double total_weight;
    pf_sample_t *sample;
    pf_vector_t pose, upper_pose;
    pf_vector_t hit;

    bool isColorSeen = false;


    self = (AMCLLaser*) data->sensor;

    total_weight = 0.0;

    /** Following a model similar to Likelihood Field when we have a colored blob sighting */
    // Pre-compute a couple of things
    double z_hit_denom = 2 * self->sigma_hit * self->sigma_hit;
    double z_rand_mult = 1.0/data->range_max;


    // Compute the sample weights
    for (j = 0; j < set->sample_count; j++)
    {
        sample = set->samples + j;
        pose = sample->pose;

        // Take account of the laser pose relative to the robot
        upper_pose = pf_vector_coord_add(self->laser_pose, upper_sample->pose);

        p = 1.0;
        pz = 0.0;
        pz_color = 0.0;

        step = (data->range_count - 1) / (self->max_beams - 1);
        for (i = 0; i < data->range_count; i += step)
        {
            obs_range = data->ranges[i][0];
            obs_bearing = data->ranges[i][1];
            obs_color = data->ranges[i][2];

            int x0, y0;
            x0 = MAP_GXWX(self->map, pose.v[0]);
            y0 = MAP_GYWY(self->map, pose.v[1]);

            if( !MAP_VALID(self->map, x0, y0) || (self->map->cells[MAP_INDEX(self->map,x0,y0)].occ_state > -1) ){
                z = self->map->max_occ_dist;
                //color_z = self->color_map->max_occ_dist;
            }

            else{
                // Compute the endpoint of the beam
                hit.v[0] = upper_pose.v[0] + obs_range * cos(upper_pose.v[2] + obs_bearing);
                hit.v[1] = upper_pose.v[1] + obs_range * sin(upper_pose.v[2] + obs_bearing);

                if(obs_color == 50){

                    //ROS_INFO("inside colored weighting");

                    int x1, y1;

                    x1 = MAP_GXWX(self->map, hit.v[0]);
                    y1 = MAP_GYWY(self->map, hit.v[1]);



                    // Part 1: Get distance from the hit to sample.
                    z = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)) * self->map->scale;

                    ROS_INFO("obs_range: %0.3f obs_bearing: %0.3f ---- sample_pose | %0.3f, %0.3f --- upper_pose | %0.3f, %0.3f --- hit_pose | %0.3f, %0.3f --- z: %0.3f"
                             ,obs_range, obs_bearing, pose.v[0], pose.v[1], upper_pose.v[0], upper_pose.v[1], hit.v[0], hit.v[1], z);

                    isColorSeen = true;
                }
                else{
                    z = self->map->max_occ_dist - 0.1;
                }
            }

            if(obs_color != 50){
                // Gaussian model
                // NOTE: this should have a normalization of 1/(sqrt(2pi)*sigma)
                pz = pz + (self->z_hit * exp(-(z * z) / z_hit_denom));

                // Part 2: random measurements
                pz = pz + (self->z_rand * z_rand_mult) *(self->z_rand * z_rand_mult);
            }
            else{
                // Gaussian model
                // NOTE: this should have a normalization of 1/(sqrt(2pi)*sigma)
                pz_color = pz_color + (self->z_hit * exp(-(z * z) / z_hit_denom));

                // Part 2: random measurements
                pz_color = pz_color + (self->z_rand * z_rand_mult) *(self->z_rand * z_rand_mult);
            }

        }
    }
    //    if(pz > 1.0){
    //        pz = 1.0;
    //    }
    //    else if(pz<0.0){
    //        pz = 0.0;
    //    }

    assert(pz <= 1.0);
    assert(pz >= 0.0);

    assert(pz_color <= 1.0);
    assert(pz_color >= 0.0);

    if(isColorSeen){
        //      p *= pz;
        // here we have an ad-hoc weighting scheme for combining beam probs
        // works well, though...
        p += pz_color*pz_color*pz_color;
    }
    else{
        //      p *= pz;
        // here we have an ad-hoc weighting scheme for combining beam probs
        // works well, though...
        p += pz*pz*pz;
    }

    sample->weight *= p;
    total_weight += sample->weight;


    return(total_weight);
}
