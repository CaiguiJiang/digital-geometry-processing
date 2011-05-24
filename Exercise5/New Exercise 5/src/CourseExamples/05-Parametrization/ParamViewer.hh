//=============================================================================
//
//  CLASS ParamViewer
//
//=============================================================================


#ifndef PARAMETRIZATION_VIEWER_HH
#define PARAMETRIZATION_VIEWER_HH


//== INCLUDES =================================================================

#include "MeshViewer.hh"
#include <gmm.h>

//== CLASS DEFINITION =========================================================

	      

class ParamViewer : public MeshViewer
{
public:
  enum ParameterizationMode {NoParameterization, Uniform, Harmonic}; 

  typedef gmm::dense_matrix<double>  gmmMatrix;
	typedef std::vector<double>        gmmVector;

   
  /// default constructor
  ParamViewer(const char* _title, int _width, int _height, int iRepeat=1);

  // destructor
  ~ParamViewer();

  /// open mesh
  virtual bool open_mesh(const char* _meshfilename);


private:
  void calc_harmonic_parameterization();

  void calc_harmonic_parameterization_aux();

  void calc_uniform_parameterization_aux(gmmVector& x, gmmVector& y);
  
  void calc_uniform_parameterization();

  void calc_distortion(ParameterizationMode imode);

  void ComputeTextureCoordinates(int iTextureWidth, int iTextureHeight, int iReapeats, ParameterizationMode imode);

//  void calc_uniform_boundary_weights (Mesh::HalfedgeHandle heh_cur, Mesh::HalfedgeHandleheh_boundary, int boundry_len);

//  void calc_uniform_non_boundary_weights();



  // solve linear system _A * _x = _b
	void solve_linear_system( gmmMatrix& _A, 
														gmmVector& _b, 
														gmmVector& _x );

  virtual void init();
  virtual void draw(const std::string& _draw_mode);

  virtual void keyboard(int key, int x, int y);

private:
  int _Repeat;

  OpenMesh::VPropHandleT<OpenMesh::Vec2f>   vparam_u_, vparam_h_, texcoord_u_, texcoord_h_ ;
  OpenMesh::VPropHandleT<Mesh::Point>       vpos_;
  OpenMesh::VPropHandleT<Mesh::Scalar>      vparam_index_;
  OpenMesh::EPropHandleT<Mesh::Scalar>      eweight_;

  GLuint  textureID_;
  float *_TextureCoordinates_u, *_TextureCoordinates_h;

  Mesh::Point _bbMin3D, _bbMax3D, _bbMin2D, _bbMax2D;

  ParameterizationMode _ParameterizationMode_;
};


//=============================================================================
#endif // HARMONIC_MAP_VIEWER_HH defined
//=============================================================================

