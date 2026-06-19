#include "Hessian.h"

#include <cstdio>
#include <stdexcept>
#include <ostream>

Hessian::Hessian(const std::string& filename) : filename_(filename) {
    FILE* fp = std::fopen(filename.c_str(), "r");
    if (!fp) {
        throw std::runtime_error("Cannot open Hessian file: " + filename);
    }

    if (std::fscanf(fp, "%d", &natom_) != 1 || natom_ <= 0) {
        std::fclose(fp);
        throw std::runtime_error("Could not read natom from Hessian file: " + filename);
    }

    ncoord_ = 3 * natom_;
    data_.reserve(static_cast<std::size_t>(ncoord_) * static_cast<std::size_t>(ncoord_));

    double x = 0.0;
    while (std::fscanf(fp, "%lf", &x) == 1) {
        data_.push_back(x);
    }
    std::fclose(fp);

    const std::size_t expected = static_cast<std::size_t>(ncoord_) * static_cast<std::size_t>(ncoord_);
    if (data_.size() < expected) {
        throw std::runtime_error("Not enough Hessian values in file: " + filename);
    }
    data_.resize(expected);
}

int Hessian::natom() const { return natom_; }
int Hessian::ncoord() const { return ncoord_; }

double Hessian::at(int i, int j) const {
    return data_.at(static_cast<std::size_t>(i) * static_cast<std::size_t>(ncoord_) + static_cast<std::size_t>(j));
}

void Hessian::printSummary(std::ostream& os) const {
    os << "Hessian(" << filename_ << "), natom=" << natom_ << ", ncoord=" << ncoord_;
}

std::string Hessian::className() const { return "Hessian"; }
