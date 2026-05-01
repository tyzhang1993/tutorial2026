#ifndef DATAOBJECT_H
#define DATAOBJECT_H

#include <iosfwd>
#include <string>

class DataObject {
public:
    virtual ~DataObject() = default;
    virtual void printSummary(std::ostream& os) const = 0;
    virtual std::string className() const = 0;
};

#endif
