#include "RBFEvaluator.h"

	RBFEvaluator::RBFEvaluator()
	{
		rbf = new TriHarmonic();
		type = TriHarmonic_t;
		Beta = 0.61;
	}

	RBFEvaluator::~RBFEvaluator()
	{
		delete rbf;
	}

	void RBFEvaluator::Next()
	{	
		type = (RBFType)((int)type + 1);
		if (type > B_Spline_t) type = (RBFType)0;
	}

	const char* RBFEvaluator::ToString()
	{
		switch (type)
		{
			case TriHarmonic_t: return "Triharmonic";
			case B_Spline_t: return "BSpline";
			default:
				throw std::exception("Unknown type: " + type);
		}
	}

	double RBFEvaluator::Evaluate(double r)
	{
		delete rbf;
		switch (type)
		{
			case TriHarmonic_t:	
				rbf = new TriHarmonic();
				break;
			case B_Spline_t:
				rbf = new B_Spline(Beta);
				break;
			default:
				throw std::exception("Unknown type: " + type);
		}
		return rbf->Evaluate(r);
	}
