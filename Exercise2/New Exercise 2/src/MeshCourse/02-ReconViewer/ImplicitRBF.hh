//=============================================================================


#ifndef RBF_HH
#define RBF_HH


//=============================================================================


#include <OpenMesh/Core/Math/VectorT.hh>
#include <vector>
#include <gmm.h>
#include "Implicit.h"
#include "RBFEvaluator.h"

//=============================================================================


class ImplicitRBF: public Implicit
{
public:
	
	typedef OpenMesh::Vec3f            Vec3f;
	typedef OpenMesh::Vec3d            Vec3d;
    typedef gmm::dense_matrix<double>  gmmMatrix;
	typedef std::vector<double>        gmmVector;
	
  // fit RBF to given constraints
  ImplicitRBF( const std::vector<Vec3f>& _points, 
               const std::vector<Vec3f>& _normals,
               float& epsilon,
			   RBFEvaluator* rbfEval);


  // evaluate RBF at position _p
  double operator()(const Vec3f& _p) const;

  

private:
	RBFEvaluator* rbfEvaluator;

	// evaluate basis function of RBF
	double kernel(const Vec3d& _c, const Vec3d& _x) const
  {
    double r = (_x-_c).norm();
	return rbfEvaluator->Evaluate(r);
  }  

	// solve linear system _A * _x = _b
	void solve_linear_system( gmmMatrix& _A, 
                            gmmVector& _b, 
                            gmmVector& _x );    

private:
  std::vector<Vec3d>   centers_;
  std::vector<double>  weights_;
};


//=============================================================================
#endif // RBF_HH defined
//=============================================================================

