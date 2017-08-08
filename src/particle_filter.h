/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

struct Particle {

	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};


class ParticleFilter {

	// Number of particles to draw
	int num_particles;



	// Flag, if filter is initialized
	bool is_initialized;

	// Vector of weights of all particles
	std::vector<double> weights;

public:

	// Set of current particles
	std::vector<Particle> particles;

	// Constructor
	// @param M Number of particles
	ParticleFilter() : num_particles(0), is_initialized(false) {}

	// Destructor
	~ParticleFilter() {}

	/**
	 * Initialize particle filter by initializing particles to Gaussian distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
	void init(double x, double y, double theta, double std[]);

	/**
	 * Predict the state for the next time step using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);

	/**
	 * Updates the weights for each particle based on the likelihood of the observed measurements.
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [standard deviation of range [m],
	 *   standard deviation of bearing [rad]]
	 * @param observations Vector of landmark observations
	 * @param map_landmarks Map class containing map landmarks
	 */
	void updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations,
			Map map_landmarks);

	/**
	 * Resample from the updated set of particles to form the new set of particles.
	 */
	void resample();

	/*
	 * Set a particle's list of associations, along with the associations calculated world x, y coordinates.
   * @param particle The particle whose associations are to be set
   * @param associations The landmark ID that goes along with each listed association
 	 * @param sense_x The associations x mapping already converted to world coordinates
 	 * @param sense_y The associations y mapping already converted to world coordinates
	 */
	Particle SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y);

	std::string getAssociations(Particle best);
	std::string getSenseX(Particle best);
	std::string getSenseY(Particle best);

	/**
	 * Return whether particle filter is initialized yet or not.
	 */
	const bool initialized() const {
		return is_initialized;
	}
};

#endif /* PARTICLE_FILTER_H_ */
