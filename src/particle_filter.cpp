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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 20;

	std::default_random_engine gen;

	std::normal_distribution<double> N_x(x, std[0]);
	std::normal_distribution<double> N_y(y, std[1]);
	std::normal_distribution<double> N_theta(theta, std[2]);

	for(int i=0; i<num_particles; i++){
		// create a particle
		Particle particle;
		// assign values to the particle
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;
		// append the particle to the list of particles
		particles.push_back(particle);
		weights.push_back(1);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	std::default_random_engine gen;

	for (int i = 0; i < num_particles; i++)	{
		double new_x;
		double new_y;
		double new_theta;

		if (yaw_rate == 0)
		{
			new_x = particles[i].x + velocity*cos(particles[i].theta)*delta_t;
			new_y = particles[i].y + velocity*sin(particles[i].theta)*delta_t;
			new_theta = particles[i].theta;
		}
		else {
			new_x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			new_y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			new_theta = particles[i].theta + yaw_rate * delta_t;
		}

		std::normal_distribution<double> N_x(new_x, std_pos[0]);
		std::normal_distribution<double> N_y(new_y, std_pos[1]);
		std::normal_distribution<double> N_theta(new_theta, std_pos[2]);

		particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);
		particles[i].weight = 1.0;
	}
	// cout << "weight after prediction " << particles[0].weight << particles[1].weight << particles[2].weight << endl;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

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

	std::vector<LandmarkObs> transformed_obs;	// vector to hold the transformed observations

	double distance;
	double previous_dist = 9999.99;
	
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;

	double factor = 1.0/(2.0 * M_PI * std_landmark[0] * std_landmark[1]);
	// cout << factor;
	double param;
	double sig_x_2 = 2.0 * std_landmark[0] * std_landmark[0];
	double sig_y_2 = 2.0 * std_landmark[1] * std_landmark[1];

	for (unsigned int i = 0; i < particles.size(); ++i)	// loop over all the particles
	{
		transformed_obs.clear();
		associations.clear();
		sense_x.clear();
		sense_y.clear();
		particles[i].weight = 1.0;

		for (unsigned int j = 0; j < observations.size(); ++j)	// loop over number of observations
		{
			// transform the observations to map coordinates wrt current particle
			transformed_obs.push_back(LandmarkObs());
			transformed_obs[j].x = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
			transformed_obs[j].y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);

			// build associations
			
			associations.push_back(0);	// this needes to be set to the id of the closest landmark
			// sense_x.push_back(transformed_obs[j].x);
			// sense_y.push_back(transformed_obs[j].y);

			// loop nover all the landmark & get the id of nearest landmark 
			previous_dist = 9999.99;
			for (unsigned int l = 0; l < map_landmarks.landmark_list.size(); ++l)
			{
				// calculate the distance with landmark
				distance = sqrt( (map_landmarks.landmark_list[l].x_f - transformed_obs[j].x) * 
								 (map_landmarks.landmark_list[l].x_f - transformed_obs[j].x) + 
							     (map_landmarks.landmark_list[l].y_f - transformed_obs[j].y) * 
							     (map_landmarks.landmark_list[l].y_f - transformed_obs[j].y) );
				if (distance < previous_dist)
				{
					previous_dist = distance;
					associations[j] = l;	//map_landmarks[l].id;
				}
			}
			// cout << transformed_obs[j].y - map_landmarks.landmark_list[associations[j]].y_f << " ";

			//calculate weight for the particle
			param = ( (transformed_obs[j].x - map_landmarks.landmark_list[associations[j]].x_f) *
			          (transformed_obs[j].x - map_landmarks.landmark_list[associations[j]].x_f) / sig_x_2 ) + 
					( (transformed_obs[j].y - map_landmarks.landmark_list[associations[j]].y_f) *
					  (transformed_obs[j].y - map_landmarks.landmark_list[associations[j]].y_f) / sig_y_2 );
			// cout << param << "  ";
			particles[i].weight *= factor * exp(-param);

		}
		// cout << "weight of P" << i << " = " << particles[i].weight << endl;
		weights[i] = particles[i].weight;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<Particle> new_particles;

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> distribution(weights.begin(), weights.end());

    for(int i = 0; i < num_particles; i++){
        Particle p = particles[distribution(gen)];
        new_particles.push_back(p);
    }
    particles = new_particles;

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
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
