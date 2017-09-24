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
    num_particles = 100;
    particles.resize(num_particles);
    weights.resize(num_particles, 1);

    float std_x = std[0];
    float std_y = std[1];
    float std_theta = std[2];

    // Init gaussian distributions for x, y, theta
    default_random_engine gen;

    // std::normal_distribution<type> var(dist_mean, std_dev);
    normal_distribution<double> dist_x(x, std_x);
    normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_theta(theta, std_theta);

    // initialize all particles
    for(int i = 0; i < num_particles; i++){
        float init_x = dist_x(gen);
        float init_y = dist_y(gen);
        float init_theta = dist_theta(gen);

        Particle p;

        p.id = i;
        p.x = init_x;
        p.y = init_y;
        p.theta = init_theta;
        p.weight = weights[i];

        particles[i] = p;
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    /* For yaw_rate != 0:
     * x_new = x + (velocity/yaw_rate) * (sin(theta + yaw_rate * dt) - sin(theta))
     * y_new = y + (velocity/yaw_rate) * (cos(theta) - cos(theta + yaw_rate * dt))
     * theta_new = theta + yaw_rate * dt
     * */

    float std_x = std_pos[0];
    float std_y = std_pos[1];
    float std_theta = std_pos[2];

    default_random_engine gen;

    for(int i = 0; i < num_particles; i++){
        Particle p = particles[i];

        // add gaussian noise
        normal_distribution<double> dist_x(p.x, std_x);
        normal_distribution<double> dist_y(p.y, std_y);
        normal_distribution<double> dist_theta(p.theta, std_theta);

        float noise_x = dist_x(gen);
        float noise_y = dist_y(gen);
        float noise_theta = dist_theta(gen);

        // Predict location
        float x_new;
        float y_new;
        float theta_new;

        if (abs(yaw_rate) < 0.0000001){
            x_new = noise_x + velocity * delta_t * cos(theta_new);
            y_new = noise_y + velocity * delta_t * sin(theta_new);
            theta_new = noise_theta;
        } else {
            x_new = noise_x + (velocity / yaw_rate) * (sin(noise_theta + (yaw_rate * delta_t)) - sin(noise_theta));
            y_new = noise_y + (velocity / yaw_rate) * (cos(noise_theta) - cos(noise_theta + (yaw_rate * delta_t)));
            theta_new = noise_theta + (yaw_rate * delta_t);
        }

        Particle p_new;

        p_new.x = x_new;
        p_new.y = y_new;
        p_new.theta = theta_new;
        p_new.associations = p.associations;
        p_new.sense_x = p.sense_x;
        p_new.sense_y = p.sense_y;

        particles[i] = p_new;
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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

    /* Multivariate normal distribution function: updated_weights = (exp(-1/2(x - mu).transpose() * sigma^-1(x-mu)) / sqrt(abs(2*M_PI*sigma))
     * x = i'th landmark measurement for one particular particle
     * mu = predicted measurement for the map landmark corresponding to the ith measurement
     * sigma = [sig_xx, 0,
     *          0     , sig_yy];
     *
     * */

    float std_x = std_landmark[0];
    float std_y = std_landmark[1];
    float denominator = 2 * M_PI * std_x * std_y;

    for(int i = 0; i < num_particles; i++){
        Particle p = particles[i];

        float weight_new = 1.;
        for(int o = 0; o < observations.size(); o++){
            LandmarkObs obs = observations[o];

            // equations for transforming from vehicle coordinates to map coordinates
            // xm = p.x + (cos(p.theta) * obs.x - (sin(p.theta)*obs.y);
            // ym = p.y + (sin(p.theta) * obs.x)+ (cos(p.theta)*obs.y);
            float xm = p.x + (cos(p.theta) * obs.x - (sin(p.theta) * obs.y));
            float ym = p.y + (sin(p.theta) * obs.x)+ (cos(p.theta) * obs.y);

            // pair nearest landmark with observation
            float min_dist = sensor_range;
            float x_delta = sensor_range;
            float y_delta = sensor_range;

            for(int L = 0; L < map_landmarks.landmark_list.size(); L++){
                Map::single_landmark_s landmark = map_landmarks.landmark_list[L];

                float dist_to_landmark = dist(xm,ym, landmark.x_f, landmark.y_f);

                if(dist_to_landmark < min_dist){
                    min_dist = dist_to_landmark;
                    x_delta = xm - landmark.x_f;
                    y_delta = ym - landmark.y_f;
                }
            }
            // P(x,y) = 1/(2*pi*std_x * std_y) * exp(-(x-mu_x)^2/(2*std_x^2) + (y + mu_y)^2 / (2*std_y^2))
            //
            float exponent = -((pow(x_delta, 2) / (2 * std_x * std_x)) + (pow(y_delta, 2) / (2 * std_y * std_y)));
            weight_new *= exp(exponent) / denominator;

        }
        weights[i] = weight_new;
        particles[i].weight = weight_new;

    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	// http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    // create vector for new particles of size num_particles
    vector<Particle> p_new(num_particles);

    // using example from discrete_distribution reference...
    random_device rd;
    mt19937 gen(rd());
    // Fill new particle vector to the same size as previous particle vector.
    // Particles selected proportionally their weight in the weight array.
    for(int n = 0; n < num_particles; n++){
        // discrete_distribution::(constructor) => distribution(iter.first(), iter.last())
        discrete_distribution<> distribution(weights.begin(), weights.end());
        int rand_particle = distribution(gen);
        // Fill in p_new with one of the particles from the particles[] vector.
        p_new[n] = particles[rand_particle];
    }
    // Overwite particle vector with resampled vector
    particles = p_new;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	// particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
