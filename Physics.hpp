#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include <cmath>
#include <vector>
#include "emp/math/Random.hpp"
#include "Particle.hpp"

namespace Physics {

inline double WrapDiff(double val, double max_val) {
        if (val > max_val / 2.0) return val - max_val;
        if (val < -max_val / 2.0) return val + max_val;
        return val;
    }

inline void ApplyLinearSpringForces(std::vector<Particle>& particles, double canvas_width, double canvas_height, double spring_k) {
    for (auto& p1 : particles) {
        for (size_t bonded_id : p1.GetBondedIDs()) {
            if (bonded_id >= particles.size() || p1.GetID() >= bonded_id) continue;
            auto& p2 = particles[bonded_id];

            double dx = WrapDiff(p2.GetX() - p1.GetX(), canvas_width);
            double dy = WrapDiff(p2.GetY() - p1.GetY(), canvas_height);

            double dist = std::hypot(dx, dy);
            double target_dist = p1.GetRadius() + p2.GetRadius();

            if (dist > 0.001) {
                double force = (dist - target_dist) * spring_k;
                double fx = (dx / dist) * force;
                double fy = (dy / dist) * force;

                p1.Move(fx, fy, canvas_width, canvas_height);
                p2.Move(-fx, -fy, canvas_width, canvas_height);
            }
        }
    }
}

inline void ApplyAngularSpringForces(std::vector<Particle>& particles, double canvas_width, double canvas_height, double angular_k = 0.03, double theta_0 = 140.0 * M_PI / 180.0) {
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p2 = particles[i];
        if (p2.GetSpecies() != Species::Bonding) continue;

        const auto& bonds = p2.GetBondedIDs();

        if (bonds.size() == 2) {
            size_t id1 = bonds[0];
            size_t id3 = bonds[1];
            if (id1 >= particles.size() || id3 >= particles.size()) continue;

            auto& p1 = particles[id1];
            auto& p3 = particles[id3];

            double dx1 = WrapDiff(p1.GetX() - p2.GetX(), canvas_width);
            double dy1 = WrapDiff(p1.GetY() - p2.GetY(), canvas_height);
            double dx3 = WrapDiff(p3.GetX() - p2.GetX(), canvas_width);
            double dy3 = WrapDiff(p3.GetY() - p2.GetY(), canvas_height);

            double len1 = std::hypot(dx1, dy1);
            double len3 = std::hypot(dx3, dy3);
            if (len1 < 0.001 || len3 < 0.001) continue;

            double u1x = dx1 / len1, u1y = dy1 / len1;
            double u3x = dx3 / len3, u3y = dy3 / len3;

            double a1 = std::atan2(dy1, dx1);
            double a3 = std::atan2(dy3, dx3);
            double delta_a = a3 - a1;

            while (delta_a > M_PI) delta_a -= 2.0 * M_PI;
            while (delta_a <= -M_PI) delta_a += 2.0 * M_PI;

            double current_theta = std::abs(delta_a);
            double angular_error = theta_0 - current_theta;
            double force_mag = angular_k * angular_error;

            double sign = (delta_a >= 0.0) ? 1.0 : -1.0;

            double p1x = -u1y, p1y = u1x;
            double p3x = -u3y, p3y = u3x;

            double fx1 = -sign * force_mag * p1x;
            double fy1 = -sign * force_mag * p1y;
            double fx3 =  sign * force_mag * p3x;
            double fy3 =  sign * force_mag * p3y;

            p1.Move(fx1, fy1, canvas_width, canvas_height);
            p3.Move(fx3, fy3, canvas_width, canvas_height);
            p2.Move(-(fx1 + fx3), -(fy1 + fy3), canvas_width, canvas_height);
        }
    }
}

inline void ApplyBondingAttraction(std::vector<Particle>& particles, double canvas_width, double canvas_height, double attraction_k = 0.015, double template_k = 0.035, double cutoff_dist = 14.0) {
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Species s1 = particles[i].GetSpecies();
            Species s2 = particles[j].GetSpecies();

            if (s1 == Species::Bonding && s2 == Species::Bonding) {
                if (particles[i].GetBondedIDs().size() >= 2) continue;
                if (particles[j].GetBondedIDs().size() >= 2) continue;
                if (particles[i].IsBondedWith(j)) continue;

                double dx = WrapDiff(particles[j].GetX() - particles[i].GetX(), canvas_width);
                double dy = WrapDiff(particles[j].GetY() - particles[i].GetY(), canvas_height);
                double dist = std::hypot(dx, dy);
                double min_dist = particles[i].GetRadius() + particles[j].GetRadius();

                if (dist > min_dist && dist < cutoff_dist) {
                    double force = attraction_k * (1.0 - (dist - min_dist) / (cutoff_dist - min_dist));
                    double fx = (dx / dist) * force;
                    double fy = (dy / dist) * force;
                    particles[i].Move(fx, fy, canvas_width, canvas_height);
                    particles[j].Move(-fx, -fy, canvas_width, canvas_height);
                }
            }
            else if ((s1 == Species::Catalyst && s2 == Species::Bonding) ||
                     (s1 == Species::Bonding && s2 == Species::Catalyst)) {
                
                double dx = WrapDiff(particles[j].GetX() - particles[i].GetX(), canvas_width);
                double dy = WrapDiff(particles[j].GetY() - particles[i].GetY(), canvas_height);
                double dist = std::hypot(dx, dy);
                double min_dist = particles[i].GetRadius() + particles[j].GetRadius();
                double cat_cutoff = min_dist + 10.0;

                if (dist > min_dist && dist < cat_cutoff) {
                    double force = template_k * (1.0 - (dist - min_dist) / (cat_cutoff - min_dist));
                    double fx = (dx / dist) * force;
                    double fy = (dy / dist) * force;
                    particles[i].Move(fx, fy, canvas_width, canvas_height);
                    particles[j].Move(-fx, -fy, canvas_width, canvas_height);
                }
            }
        }
    }
}

inline void ResolveRepulsion(std::vector<Particle>& particles, double canvas_width, double canvas_height, double repulsion_k) {
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            if (particles[i].IsColliding(particles[j])) {
                double dx = WrapDiff(particles[j].GetX() - particles[i].GetX(), canvas_width);
                double dy = WrapDiff(particles[j].GetY() - particles[i].GetY(), canvas_height);
                double dist = std::hypot(dx, dy);
                double min_dist = particles[i].GetRadius() + particles[j].GetRadius();

                if (dist > 0.001 && dist < min_dist) {
                    double overlap = min_dist - dist;
                    double fx = (dx / dist) * overlap * repulsion_k;
                    double fy = (dy / dist) * overlap * repulsion_k;

                    particles[i].Move(-fx, -fy, canvas_width, canvas_height);
                    particles[j].Move(fx, fy, canvas_width, canvas_height);
                }
            }
        }
    }
}

inline void ThermalBondBreaking(std::vector<Particle>& particles, double temperature, double bond_break_prob, emp::Random& random) {
    double p_break = bond_break_prob * temperature;
    for (auto& p1 : particles) {
        auto bonded_ids = p1.GetBondedIDs();
        for (size_t bonded_id : bonded_ids) {
            if (p1.GetID() < bonded_id && random.P(p_break)) {
                p1.RemoveBond(bonded_id);
                particles[bonded_id].RemoveBond(p1.GetID());
            }
        }
    }
}

} // namespace Physics

#endif