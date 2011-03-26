//=============================================================================


#include "ImplicitMLS.hh"
#include <float.h>
#include <iostream>

//== IMPLEMENTATION ==========================================================
ImplicitMLS::ImplicitMLS( const std::vector<Vec3f>& _points, 
                          const std::vector<Vec3f>& _normals )
: points_(_points), normals_(_normals)
{
  double beta = CalculateBeta();
  InvBetaSquare_ = 1.0 / pow(beta, 2.0);
}


double ImplicitMLS::operator()(const Vec3f& _p) const
{
  double dist(0);
  int i = 0;
  for (std::vector<Vec3f>::const_iterator it = points_.begin(); it != points_.end(); it++, i++)
  {
		Vec3f currentPoint = *it;
		if (_p == currentPoint) continue;
		Vec3f diffVector = (_p - currentPoint);
		double dotProduct = normals_[i]|diffVector;
		double distance = diffVector.norm();
		dist += (dotProduct * Phi(distance));
  }
  return dist;
}

double ImplicitMLS::CalculateBeta()
{
	// Calculating the average distance of each point to its closest neighbor
	//std::fstream out = std::fstream("D:\\Log.txt", std::ios_base::
	double sum = 0.0;
	for (std::vector<Vec3f>::const_iterator it = points_.begin(); it != points_.end(); it++)
	{
		Vec3f currentPoint = *it;
		double d = CalculateDistanceToClosestPoint(currentPoint);
		sum += d;
	}

	sum /= points_.size();
	return 2 * sum;
}

double ImplicitMLS::CalculateDistanceToClosestPoint(const Vec3f& p)
{
	double smallestDistance = DBL_MAX;
	for (std::vector<Vec3f>::const_iterator it = points_.begin(); it != points_.end(); it++)
	{
		Vec3f currentPoint = *it;
		if (p == currentPoint) continue;
		double distance = (p - currentPoint).norm();
		smallestDistance = (smallestDistance < distance) ? smallestDistance : distance;
	}
	return smallestDistance;
}

double ImplicitMLS::Phi(double r) const
{
	return exp(-pow(r,2.0) * InvBetaSquare_);
}
//=============================================================================
