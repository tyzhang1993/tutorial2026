#ifndef HESSIAN_H
#define HESSIAN_H

#include "DataObject.h"
#include <string>
#include <vector>

class Hessian : public DataObject {
public:
    explicit Hessian(const std::string& filename);

    int natom() const;
    int ncoord() const;
    double at(int i, int j) const;

    void printSummary(std::ostream& os) const override;
    std::string className() const override;

private:
    std::string filename_;
    int natom_ = 0;
    int ncoord_ = 0;
    std::vector<double> data_;
};

#endif
