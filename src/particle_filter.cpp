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
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	// x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.

  num_particles = 100;

  std::default_random_engine gen;

  std::normal_distribution<double> N_x(x, std[0]);
  std::normal_distribution<double> N_y(y, std[1]);
  std::normal_distribution<double> N_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++) {
    Particle particle;
    particle.id = i;
    particle.x = N_x(gen);
    particle.y = N_y(gen);
    particle.theta = N_theta(gen);
    particle.weight = 1;

    particles.push_back(particle);
    weights.push_back(1);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
  default_random_engine gen;

  for (int i = 0; i < num_particles; i++) {
    double new_x;
    double new_y;
    double new_theta;

    if (yaw_rate == 0) {
      new_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      new_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
      new_theta = particles[i].theta;
    } else {
      new_x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
      new_y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      new_theta = particles[i].theta + yaw_rate * delta_t;
    }

    normal_distribution<double> N_x(new_x, std_pos[0]);
    normal_distribution<double> N_y(new_y, std_pos[1]);
    normal_distribution<double> N_theta(new_theta, std_pos[2]);

    particles[i].x = N_x(gen);
    particles[i].y = N_y(gen);
    particles[i].theta = N_theta(gen);
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian distribution.
  for (int p = 0; p < num_particles; p++) {
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;

    vector<LandmarkObs> trans_observations;
    LandmarkObs obs;

    // Transform all observations from vehicle into map coordinates (rotation and translation)
    for (int i = 0; i < observations.size(); i++) {
      LandmarkObs trans_obs;
      obs = observations[i];

      trans_obs.x = particles[p].x + (obs.x * cos(particles[p].theta) - obs.y * sin(particles[p].theta));
      trans_obs.y = particles[p].y + (obs.x * sin(particles[p].theta) + obs.y * cos(particles[p].theta));
      trans_observations.push_back(trans_obs);
    }

    particles[p].weight = 1.0;

    // Associate observations with landmarks based on Euclidean distance
    for (int i = 0; i < trans_observations.size(); i++) {

      double closest_dist = sensor_range;
      int association = 0;

      for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
        double landmark_x = map_landmarks.landmark_list[j].x_f;
        double landmark_y = map_landmarks.landmark_list[j].y_f;

        double distance = dist(landmark_x, landmark_y, trans_observations[i].x, trans_observations[i].y);

        if (distance < closest_dist) {
          closest_dist = distance;
          association = j;
        }
      }

      if(association != 0) {
        double meas_x = trans_observations[i].x;
        double meas_y = trans_observations[i].y;
        double mu_x = map_landmarks.landmark_list[association].x_f;
        double mu_y = map_landmarks.landmark_list[association].y_f;

        long double factor = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]) * exp(-((pow(meas_x - mu_x, 2) / (2 * std_landmark[0] * std_landmark[0]) + pow(meas_y - mu_y, 2) / (2 * std_landmark[1] * std_landmark[1]))));

        // Product
        if (factor > 0) {
          particles[p].weight *= factor;
        }
      }

      associations.push_back(association + 1);
      sense_x.push_back(trans_observations[i].x);
      sense_y.push_back(trans_observations[i].y);
    }

    particles[p] = SetAssociations(particles[p], associations, sense_x, sense_y);
    weights[p] = particles[p].weight;
  }
}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight.

  default_random_engine gen;
  discrete_distribution<int> distribution(weights.begin(), weights.end());

  vector<Particle> resample_particles;

  for (int i = 0; i < num_particles; i++) {
    resample_particles.push_back(particles[distribution(gen)]);
  }

  particles = resample_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y) {
	// Set a particles list of associations, along with the associations calculated world x, y coordinates.

	// Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations = associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best) {
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // Get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best) {
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // Get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best) {
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // Get rid of the trailing space
    return s;
}
