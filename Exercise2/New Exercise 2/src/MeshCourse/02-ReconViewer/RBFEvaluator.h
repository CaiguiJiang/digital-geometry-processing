#ifndef RBFEVALUATOR_H
#define RBFEVALUATOR_H

#include <iostream>

class RBFEvaluator
{
private:
	/* Always keep B_Spline_t as the last type!!! */
	typedef enum {TriHarmonic_t, B_Spline_t} RBFType;

	class RBF
	{
		public:
			virtual double Evaluate(double r) = 0;
	};

	class TriHarmonic : public RBF
	{
		public:
			double Evaluate(double r) { return r*r*r; }
	};

	class B_Spline : public RBF
	{
		private:
			//TODO: make Beta double
			double Beta;
		public:
			B_Spline(double beta = 0.1) 
			{ 
				Beta = beta; 				
			}
			double Evaluate(double r) 
			{ 
				double s = r / Beta;
				if (s >= -2 && s < -1)
				{
					double p = 2+s;
					return p*p*p / 6;
				}
				else if (s >= -1 && s < 0)
				{
					return (4-6*s*s-3*s*s*s) / 6;
				}
				else if (s >= 0 && s < 1)
				{
					return (4-6*s*s+3*s*s*s) / 6;
				}
				else if (s >= 1 && s < 2)
				{
					double p = 2-s;
					return p*p*p / 6;
				}
				else return 0.0;
			}	
	};

	RBF* rbf;
	RBFType type;
public:
	double Beta;
	RBFEvaluator();

	~RBFEvaluator();

	//TODO: Remove {double epsilon}
	void Next();

	const char* ToString();


	double Evaluate(double r);
};

#endif