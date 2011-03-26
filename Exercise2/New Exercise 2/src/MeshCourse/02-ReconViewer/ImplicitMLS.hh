//=============================================================================


#ifndef MLS_HH
#define MLS_HH


//=============================================================================


#include <OpenMesh/Core/Math/VectorT.hh>
#include <vector>
#include <float.h>
#include "Implicit.h"


//=============================================================================


class ImplicitMLS: public Implicit
{
public:

	typedef OpenMesh::Vec3f Vec3f;

	
  // fit RBF to given constraints
  ImplicitMLS( const std::vector<Vec3f>& _points, 
               const std::vector<Vec3f>& _normals );

  // evaluate implicit at position _p
  double operator()(const Vec3f& _p) const;


private:

    const std::vector<Vec3f>&  points_;
    const std::vector<Vec3f>&  normals_;
    float                      InvBetaSquare_;

	double CalculateBeta();
	double CalculateDistanceToClosestPoint(const Vec3f& p);
	double Phi(double r) const;
};


//=============================================================================
#endif // RBF_HH defined
//=============================================================================

