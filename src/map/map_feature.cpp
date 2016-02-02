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
/**************************************************************************
 * Desc: Global map storage functions
 * Author: Andrew Howard
 * Date: 6 Feb 2003
 * CVS: $Id: map_store.c 2951 2005-08-19 00:48:20Z gerkey $
**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>

#include "map.h"

using namespace std;

int doesFileExist(const char* filename);

// SINA: Return the value for crosswalk feature based on the current pose
int map_see_crosswalk(map_t *map, pf_vector_t pose[3]){



}

////////////////////////////////////////////////////////////////////////////
// SINA: Load feature properties from the file
void map_feature_load(map_t *map, const char *filename){

    printf("Loading features from ");
    printf(filename);
    printf("\n");

    if (!doesFileExist(filename)){
        printf("The feature cannot be loaded from the provided file!\n");
        return;
    }

    // File exists! Now let's load the features...

    map->cross_walks_count = 0;
    map->turn_points_count = 0;
    map->junctions_count = 0;

    FILE* fptr = fopen(filename, "r");
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    bool read_data = false;
    while ((read = getline(&line, &len, fptr)) != -1) {

        std::string line_string(line);
        vector<string> splitted;
        stringstream ssin(line_string);
        while (ssin.good()){
            string temp;
            ssin >> temp;
            splitted.push_back(temp);
        }
        splitted.pop_back();

        if (!read_data && splitted.size() != 0 && splitted[0] == "Features:"){
            read_data = true;
            continue;
        }
        if (!read_data){
            continue;
        }

        if (splitted[0] == "CrossWalk"){
            cross_walk_t cross_walk;
            cross_walk.x = atof(splitted[1].c_str());
            cross_walk.y = atof(splitted[2].c_str());

            map->cross_walks[map->cross_walks_count] = cross_walk;
            map->cross_walks_count++;
        }else if (splitted[0] == "TurnPoint"){
            turn_point_t turn_point;
            turn_point.x = atof(splitted[1].c_str());
            turn_point.y = atof(splitted[2].c_str());

            if (splitted[3] == "LEFT")
                turn_point.orientation = 0;
            else if (splitted[3] == "RIGHT")
                turn_point.orientation = 1;

            map->turn_points[map->turn_points_count] = turn_point;
            map->turn_points_count++;
        }else if (splitted[0] == "Junction"){
            junction_t junction;
            junction.x = atof(splitted[1].c_str());
            junction.y = atof(splitted[2].c_str());

            if (splitted[3] == "THREE")
                junction.type = 3;
            else if (splitted[3] == "FOUR")
                junction.type = 4;

            map->junctions[map->junctions_count] = junction;
            map->junctions_count++;
        }
    }

//    printf("\nMap Cross Walks: \n");
//    for (size_t i = 0; i < map->cross_walks_count; i++){
//        printf("Crosswalk %d:  %f\t%f\n", (int)i, map->cross_walks[i].x, map->cross_walks[i].y);
//    }

//    printf("\nTurn Points Walks: \n");
//    for (size_t i = 0; i < map->turn_points_count; i++){
//        printf("TurnPoint %d:  %f\t%f\t%d\n", (int)i, map->turn_points[i].x, map->turn_points[i].y,
//               map->turn_points[i].orientation);
//    }

//    printf("\nJunctions: \n");
//    for (size_t i = 0; i < map->junctions_count; i++){
//        printf("Junction %d:  %f\t%f\t%d\n", (int)i, map->junctions[i].x, map->junctions[i].y,
//               map->junctions[i].type);
//    }

    fclose(fptr);
    if (line) free(line);
}

int doesFileExist(const char* filename)
{
    FILE* fptr = fopen(filename, "r");
    if (fptr != NULL)
    {
        fclose(fptr);
        return 1;
    }
    return 0;
}

