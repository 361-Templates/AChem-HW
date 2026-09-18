#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>
#include <string>
#include <algorithm>

enum class Species {
    Substrate,  
    Catalyst,   
    Bonding      
};

class Particle {
private:
    size_t id;
    double x{0.0};
    double y{0.0};
    double radius{5.0};
    double vx{0.0};
    double vy{0.0};
    Species species{Species::Substrate};
    
    std::vector<size_t> bonded_ids;
    size_t max_bonds{2}; // Maximum 2 bonds per particle allows 1D chains and rings

public:

    // Static lookup mapping each species to its physical size
    static double GetDefaultRadius(Species s) {
        switch (s) {
            case Species::Substrate: return 3.0;
            case Species::Bonding:    return 5.0;
            case Species::Catalyst:  return 12.0;
        }
        return 5.0;
    }

    Particle(size_t id, double x, double y, Species species = Species::Substrate, double radius = -1.0)
        : id(id), x(x), y(y), species(species) {
        // Default to species radius if no explicit override is passed
        this->radius = (radius > 0.0) ? radius : GetDefaultRadius(species);
    }

    // Getters
    size_t GetID() const { return id; }
    double GetX() const { return x; }
    double GetY() const { return y; }
    double GetVX() const { return vx; }
    double GetVY() const { return vy; }
    double GetRadius() const { return radius; }
    Species GetSpecies() const { return species; }
    const std::vector<size_t>& GetBondedIDs() const { return bonded_ids; }

    // Setters & Modifiers
    void SetSpecies(Species new_species) { 
        species = new_species; 
        radius = GetDefaultRadius(new_species);
    }

    void Move(double dx, double dy, double max_w, double max_h) {
        x += dx;
        y += dy;

        // Toroidal boundary wrapping
        if (x < 0.0) x += max_w;
        else if (x >= max_w) x -= max_w;

        if (y < 0.0) y += max_h;
        else if (y >= max_h) y -= max_h;
    }

    bool IsColliding(const Particle& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dist_sq = (dx * dx) + (dy * dy);
        double min_dist = radius + other.radius;
        return dist_sq <= (min_dist * min_dist);
    }

    bool CanBond() const {
        return species == Species::Bonding && bonded_ids.size() < max_bonds;
    }

    bool IsBondedWith(size_t other_id) const {
        return std::find(bonded_ids.begin(), bonded_ids.end(), other_id) != bonded_ids.end();
    }

    bool AddBond(size_t other_id) {
        if (CanBond() && !IsBondedWith(other_id)) {
            bonded_ids.push_back(other_id);
            return true;
        }
        return false;
    }

    bool RemoveBond(size_t other_id) {
        auto it = std::find(bonded_ids.begin(), bonded_ids.end(), other_id);
        if (it != bonded_ids.end()) {
            bonded_ids.erase(it);
            return true;
        }
        return false;
    }

    void ClearBonds() {
        bonded_ids.clear();
    }

    // Feel free to change these!
    std::string GetColor() const {
        switch (species) {
            case Species::Substrate: return "teal";
            case Species::Catalyst:  return "orange";
            case Species::Bonding:    return "purple";
        }
        return "black";
    }
};

#endif