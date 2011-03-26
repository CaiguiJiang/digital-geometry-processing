//=============================================================================


#include "ImplicitRBF.hh"
#include "RBFEvaluator.h"

//== IMPLEMENTATION ==========================================================


ImplicitRBF::
ImplicitRBF( const std::vector<Vec3f>& _points, 
			 const std::vector<Vec3f>& _normals,
             float& epsilon,
			 RBFEvaluator* rbfEval)
{
	rbfEvaluator = rbfEval;
	int n = _points.size();
	gmmMatrix A(2*n,2*n); // Matrix
	gmmVector b; // Right-hand-side vector (0,1)
	gmmVector x(2*n); // Variables

	for (int i = 0; i < 2*n; i++)
	{
		b.push_back((i < n) ? 0 : 1);

		for (int j = 0; j < 2*n; j++)
		{
			Vec3d x = (i < n) ? (Vec3d)_points[i] : (Vec3d)(_points[i-n] + _normals[i-n] * epsilon);
			Vec3d c = (j < n) ? (Vec3d)_points[j] : (Vec3d)(_points[j-n] + _normals[j-n] * epsilon);
			A(i,j) = kernel(x,c);
		}
	}
	
	solve_linear_system(A,b,x);
	weights_.clear();
	for (int i = 0; i < 2*n; i++)
	{
		weights_.push_back(x[i]);
		centers_.push_back(	(i < n) ? (Vec3d)_points[i] : (Vec3d)(_points[i-n] + _normals[i-n] * epsilon) );
	}
}


//-----------------------------------------------------------------------------


void
ImplicitRBF::solve_linear_system( gmmMatrix& _M, 
								  gmmVector& _b, 
								  gmmVector& _x)
{
	// solve linear system by gmm's LU factorization
	unsigned int N = _b.size();
	_x.resize(N);
  std::vector< size_t >  ipvt(N);
  gmm::lu_factor( _M, ipvt );
  gmm::lu_solve( _M, ipvt, _x, _b );
}


//-----------------------------------------------------------------------------


double 
ImplicitRBF::operator()(const Vec3f& _p) const
{

	//implement your own RBF Evaluation
  std::vector<Vec3d>::const_iterator  
    c_it(centers_.begin()),
    c_end(centers_.end());

  std::vector<double>::const_iterator   
    w_it(weights_.begin());

	const Vec3d p(_p);
  double f(0);

  for (; c_it!=c_end; ++c_it, ++w_it)
    f += *w_it * kernel(*c_it, p);

  return f;
}


//=============================================================================
