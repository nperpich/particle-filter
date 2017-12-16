/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include <stdlib.h>
#include <time.h>
#include <map> //for discrete_distribution

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 50; 

	double sig_x = std[0];
	double sig_y = std[1];
	double sig_theta = std[2];
	
	std::default_random_engine generator;
	std::normal_distribution<double> x_noisey(x, sig_x);
	std::normal_distribution<double> y_noisey(y, sig_y);
	std::normal_distribution<double> theta_noisey(theta, sig_theta);
	
	for(int i=0; i<num_particles; i++){
		//intermediate that holds all info for particle being created
		Particle sample_particle;
		
		sample_particle.id = i;
		sample_particle.x = x_noisey(generator);
		sample_particle.y = y_noisey(generator);
		sample_particle.theta = theta_noisey(generator);
		sample_particle.weight = 1.0;
		
		particles.push_back(sample_particle);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	double sig_x = std_pos[0];
	double sig_y = std_pos[1];
	double sig_theta = std_pos[2];
	
	std::default_random_engine generator;
	std::normal_distribution<double> x_noise(0.0, sig_x);
	std::normal_distribution<double> y_noise(0.0, sig_y);
	std::normal_distribution<double> theta_noise(0.0, sig_theta);
	
	for(int i=0; i<num_particles; i++){
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		
		if (abs(yaw_rate) > 0.0001) {

			particles[i].x += velocity / yaw_rate*(sin(theta + yaw_rate*delta_t) - sin(theta)) + x_noise(generator);
			particles[i].y += velocity / yaw_rate*(cos(theta) - cos(theta + yaw_rate*delta_t)) + y_noise(generator);
			particles[i].theta += yaw_rate*delta_t + theta_noise(generator);
		}
		else {
			particles[i].x += velocity*cos(theta)*delta_t + x_noise(generator);
			particles[i].y += velocity*sin(theta)*delta_t + y_noise(generator);
			particles[i].theta += theta_noise(generator);
		}
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	//double weight = 1.0;
	for (int i = 0; i< observations.size(); i++) {
		/*
		double x_obs_map = x + observations[i].x*cos(theta) - observations[i].y*sin(theta);
		double y_obs_map = y + observations[i].x*sin(theta) + observations[i].y*cos(theta);
		*/

		double lowest_dist = 160.0;
		double current_dist;
		int id = 0;
		LandmarkObs closest_landmark;
		for (int j = 0; j < predicted.size(); j++) {
			current_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if (current_dist < lowest_dist) {
				lowest_dist = current_dist;
				closest_landmark = predicted[j];
				id = j;
			
			}
		// set the observation's id to the nearest predicted landmark's id
		observations[i].id = id;
		}
	}


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	//transform observations into map coordinate system
	for(int i=0; i<num_particles; i++){
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		
		std::vector<LandmarkObs> obs_in_map;
		LandmarkObs landmarkHolder;


		for (int j = 0; j < observations.size(); j++) {
			double x_obs_map = x + observations[j].x*cos(theta) - observations[j].y*sin(theta);
			double y_obs_map = y + observations[j].x*sin(theta) + observations[j].y*cos(theta);

			landmarkHolder.id = j;
			landmarkHolder.x = x_obs_map;
			landmarkHolder.y = y_obs_map;
			
			// observations from relative to global (map cooridnates)
			obs_in_map.push_back(landmarkHolder);
		}
		std::vector<LandmarkObs> predicted;
		for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
			
			if (dist(x, y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f) <= sensor_range) {
				LandmarkObs obs1;
				obs1.id = k;
				obs1.x = map_landmarks.landmark_list[k].x_f;
				obs1.y = map_landmarks.landmark_list[k].y_f;
				predicted.push_back(obs1);
			}
		}

	    dataAssociation(predicted, obs_in_map);
		
		particles[i].weight = 1.0;
		for (int l = 0; l < obs_in_map.size(); l++) {
			double associated_x = predicted[obs_in_map[l].id].x;
			double associated_y = predicted[obs_in_map[l].id].y;
			double sig_x = std_landmark[0];
			double sig_y = std_landmark[1];
			double exp_calc = -1 * ((associated_x - obs_in_map[l].x)*(associated_x - obs_in_map[l].x) / (2 * sig_x*sig_x) + (associated_y - obs_in_map[l].y)*(associated_y - obs_in_map[l].y) / (2 * sig_y*sig_y));
			particles[i].weight *= 1 / (2 * 3.14*sig_x*sig_y)*exp(exp_calc);
		}

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	weights.clear();
	
	for (int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
	}
	std::vector<Particle> sampled_particles;

	random_device rd;
	mt19937 gen(rd());
	discrete_distribution<> d(weights.begin(), weights.end());
	map<int, int> m;
	for (int n = 0; n<num_particles; ++n) {
		sampled_particles.push_back(particles[d(gen)]);
	}
	particles = sampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<double>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<double>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
