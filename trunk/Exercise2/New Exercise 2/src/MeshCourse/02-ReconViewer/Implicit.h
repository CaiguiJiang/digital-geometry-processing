#ifndef IMPLICIT_HEADER_FILE
#define IMPLICIT_HEADER_FILE

#include <OpenMesh/Core/Math/VectorT.hh>

class Implicit{
public:
	Implicit(){}
	~Implicit(){}

	virtual double operator()(const OpenMesh::Vec3f& _p) const {return 0;}  //stub function
};


#endif