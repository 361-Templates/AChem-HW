#include "emp/web/Animate.hpp"
#include "emp/web/web.hpp"
#include "emp/math/Random.hpp"
#include "Particle.hpp"
#include "Physics.hpp"

emp::web::Document doc{"target"};

class BrownianAnimator : public emp::web::Animate {

double canvas_width = 600.0;
double canvas_height = 600.0;
size_t num_particles = 150;

// The following parameters control various aspects of the physics, feel free to tweak them and see what happens!
double spring_constant = 0.08;
double temperature = 1.8;
double bond_break_probability = 0.00005;
double spontaneous_substrate_probability = 0.01;
double repulsion_constant = 0.8;
double bonding_decay_probability = 0.0;

// where we'll draw
emp::web::Canvas canvas{canvas_width, canvas_height, "canvas"};
emp::Random random;
std::vector<Particle> particles;

public:

    BrownianAnimator() {
        // shove canvas into the div
        // along with some control buttons
        doc << canvas;
        doc << GetToggleButton("Toggle");
        doc << GetStepButton("Step");

        // Spawn particles at random initial positions
        for (size_t i = 0; i < num_particles; ++i) {
            double start_x = random.GetDouble(0.0, canvas_width);
            double start_y = random.GetDouble(0.0, canvas_height);

            Species type = Species::Substrate;
            if (i < 50) {
                type = Species::Catalyst;
            } else if (i > 80) {
                type = Species::Bonding; 
            }
            particles.emplace_back(i, start_x, start_y, type);
        }
    }

    void DegradeAndRecycle() {
        for (auto& p : particles) {
            // =========================================================
            // TODO: IMPLEMENT DECAY LOGIC HERE
            // =========================================================
        }
    }

    void ResolveCollisionsAndReactions() {
        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                // =========================================================
                // TODO: IMPLEMENT REACTION & BONDING LOGIC HERE
                // =========================================================
                
            }
        }
    }

    void DoFrame() override {
        canvas.Clear();

        if (random.P(spontaneous_substrate_probability)) {
            particles.emplace_back(0, 0, 0, Species::Substrate);
        }

        int sub_steps = 4;
        double step_temp = (temperature / std::sqrt(sub_steps)) * 0.4; // Scaled for smaller moves

        for (int step = 0; step < sub_steps; ++step) {
            // Scaled thermal movement
            for (auto& p : particles) {
                double dx = random.GetNormal(0.0, step_temp);
                double dy = random.GetNormal(0.0, step_temp);
                p.Move(dx, dy, canvas_width, canvas_height);
            }

            // Apply physics & structural forces
            Physics::ThermalBondBreaking(particles, temperature, bond_break_probability, random);
            Physics::ApplyBondingAttraction(particles, canvas_width, canvas_height);
            Physics::ApplyLinearSpringForces(particles, canvas_width, canvas_height, spring_constant);
            Physics::ApplyAngularSpringForces(particles, canvas_width, canvas_height);
            Physics::ResolveRepulsion(particles, canvas_width, canvas_height, repulsion_constant);

            // Perform chemical reactions
            DegradeAndRecycle();
            ResolveCollisionsAndReactions();
        }

        // Render bonds
        for (const auto& p1 : particles) {
            for (size_t bonded_id : p1.GetBondedIDs()) {
                if (p1.GetID() < bonded_id) {
                    const auto& p2 = particles[bonded_id];
                    canvas.Line(p1.GetX(), p1.GetY(), p2.GetX(), p2.GetY(), "purple");
                }
            }
        }

        // Render particles
        for (const auto& p : particles) {
            canvas.Circle(p.GetX(), p.GetY(), p.GetRadius(), p.GetColor(), "black");
        }


    }


};

BrownianAnimator animator;

int main() {
    //Have animator call DoFrame once to start
    animator.Step();
}