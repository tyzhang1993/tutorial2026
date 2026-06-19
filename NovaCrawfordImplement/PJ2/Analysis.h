#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <iosfwd>

class Analysis {
public:
    virtual ~Analysis() = default;
    virtual void run() = 0;
    virtual void print(std::ostream& os) const = 0;
};

#endif
