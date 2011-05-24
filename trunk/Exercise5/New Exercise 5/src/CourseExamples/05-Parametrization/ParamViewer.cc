//=============================================================================
//
//  CLASS ParamViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "StdAfx.h"
#include "ParamViewer.hh"

//== IMPLEMENTATION ========================================================== 

static bool _ParameterizationComputed_u = false, _ParameterizationComputed_h = false;
static bool _BoundingBox2DComputed = false;


ParamViewer::
ParamViewer(const char* _title, int _width, int _height, int iRepeat)
  : MeshViewer(_title, _width, _height)
{ 
  _Repeat = iRepeat*5;

  mesh_.request_vertex_colors();

  mesh_.add_property(vpos_);
  mesh_.add_property(vparam_u_);
  mesh_.add_property(vparam_h_);
  mesh_.add_property(texcoord_u_);
  mesh_.add_property(texcoord_h_);
  mesh_.add_property(eweight_);
  mesh_.add_property(vparam_index_);

  add_draw_mode("UV Domain");
	add_draw_mode("Textured mesh");

  _TextureCoordinates_u=NULL;
  _TextureCoordinates_h=NULL;

  _ParameterizationMode_ = NoParameterization;

	init();
}

ParamViewer::
~ParamViewer()
{ 
  if (glIsTexture(textureID_))  
		glDeleteTextures( 1, &textureID_);

  if (_TextureCoordinates_u)
    delete [] _TextureCoordinates_u;

  if (_TextureCoordinates_h)
    delete [] _TextureCoordinates_h;
}

void
ParamViewer::
init()
{
	// base class first
	MeshViewer::init();

  // generate checkerboard-like image
  GLubyte tex[256][256][3];
  int index=0;
  for (int x=0; x<256; ++x)
  {
    for (int y=0; y<256; ++y)
    {
      if ((x<128&&y<128) ||(x>128&&y>128))
      {
        tex[x][y][0] = 0;
        tex[x][y][1] = 255;
        tex[x][y][2] = 0;
      }
      else
      {
        tex[x][y][0] = 255;
        tex[x][y][1] = 255;
        tex[x][y][2] = 255;
      }
    }
  }
  // generate texture
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 
  if (!glIsTexture(textureID_))
    glGenTextures(1, &textureID_);
  glBindTexture(GL_TEXTURE_2D, textureID_);

  // copy texture to GL
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, 256, 256,
    0, GL_RGB, GL_UNSIGNED_BYTE, tex);
}

//-----------------------------------------------------------------------------


void
ParamViewer::solve_linear_system( gmmMatrix& _M, gmmVector& _b, gmmVector& _x)
{
  /*gmm::identity_matrix PS;   // Optional scalar product for cg
  gmm::identity_matrix PR;   // Optional preconditioner
  gmm::iteration iter(10E-8);// Iteration object with the max residu
  iter.set_maxiter(10000);
  //gmm::cg(_M, _x, _b, PS, PR, iter);
  //gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<double> > > P(_M, 10, 1e-4);
  gmm::bicgstab(_M, _x, gmm::mat_const_col(b,0), PR, iter);
  //gmm::least_squares_cg(_M, _x, _b, iter);
  //gmm::qmr(_M, _x, _b, PR, iter);
  bool ResOK = iter.converged(); */ 

  // computation of a preconditioner (ILUT)
  /*gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<double> > > P(_M, 10, 1e-4);
  gmm::iteration iter(1E-8);  // defines an iteration object, with a max residu of 1E-8
  gmm::gmres(_M, _x, _b, P, 50, iter);  // execute the GMRES algorithm
  std::cout << "The result " << _x << endl;*/ 

  unsigned int N = _b.size();
	_x.resize(N);
  std::vector< size_t >  ipvt(N);
  gmm::lu_factor( _M, ipvt );
  gmm::lu_solve( _M, ipvt, _x, _b );
} 


//-----------------------------------------------------------------------------



void ParamViewer::
calc_uniform_parameterization()
{
   // ------------- IMPLEMENT HERE ---------
	// TASK 5.1 Uniform map computation:
	// Search and order boundary vertices
  // Compute boundary parameters
  // Solve the linear system for internal vertex parameters using solve_linear_system()
  // Store the parameters in the vparam_u_ vertex property
	// ------------- IMPLEMENT HERE ---------
	for (Mesh::EdgeIter e_it=mesh_.edges_begin(); e_it!=mesh_.edges_end();++e_it)
	{
		mesh_.property(eweight_,e_it) = 2; //1 for each vetex
	}
	gmmVector a,b; 
	// Calculate the uniform parameterization.
	calc_uniform_parameterization_aux (a,b);

	//save to vertices
	for (Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it)
	{
		unsigned int i = mesh_.property (vparam_index_,v_it.handle());
		Vec2f u(a[i],b[i]);
		mesh_.property(vparam_u_,v_it) = u;
	}

}

void ParamViewer::
calc_uniform_parameterization_aux(gmmVector& x, gmmVector& y)
{
	size_t n = mesh_.n_vertices();
	gmmMatrix W(n,n); //create a matrix W in size nXn
	gmmVector bx(n); //create a vector bx in size n  
	gmmVector by(n); //create a vector by in size n
	float PI = 3.1415926535897932384626433832795;
    
	// Pass all over the edges and find an edge that on the border

    HalfedgeHandle heh_boundary;
	//bool found_boundary = false;
    for(Mesh::HalfedgeIter heh_it=mesh_.halfedges_begin(); heh_it!=mesh_.halfedges_end(); ++heh_it)
	{
		if(mesh_.is_boundary(heh_it))
		{
            heh_boundary = heh_it.handle();
			//found_boundary = true;
            break;
        }
    }

	//Now, We need to calculate what is the length of the boundary

	Mesh::HalfedgeHandle heh_cur = heh_boundary; //work with the border half edge
	Mesh::Scalar boundry_len = 0; //The length of the boundary
	int num_edges_boundary =0; // How much edges in the boundary
	do
	{
		boundry_len += mesh_.calc_edge_length(heh_cur); 
		heh_cur = mesh_.next_halfedge_handle (heh_cur);
		num_edges_boundary++;
	} while (heh_cur != heh_boundary);

	//Iterate over boundet

	heh_cur = heh_boundary;
	


	//calc_uniform_boundary_weights
	//***************************************
		
	unsigned int i = 0;
	float theta = 0;
	do
	{
		Mesh::VertexHandle vh = mesh_.to_vertex_handle (heh_cur);
		mesh_.property(vparam_index_,vh) = i;

		W(i,i) = 1; 
		bx[i] = cos(theta);
		by[i] = sin (theta);

		heh_cur = mesh_.next_halfedge_handle (heh_cur);
		i++;
		theta += (mesh_.calc_edge_length(heh_cur) / boundry_len)*2*PI;
	} while (heh_cur != heh_boundary);

	
	//find indices non-boundry vertices

	for (Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it)
	{
		if (!mesh_.is_boundary(v_it))
		{
			mesh_.property (vparam_index_,v_it.handle())=i++;
		}
	}
	
	
	//calc_uniform_non_boundary_weights
	//***************************************
	for (Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it)
	{
		if (!mesh_.is_boundary(v_it))
		{
			unsigned int i;
			i = mesh_.property(vparam_index_,v_it); //Get first vertex in Edge
			float total_weight = 0;

			//All outgoing half edges
			for (Mesh::VertexOHalfedgeIter voh_it=mesh_.voh_iter(v_it.handle()); voh_it; ++voh_it)
			{
				Mesh::VertexHandle vh_to = mesh_.to_vertex_handle (voh_it.handle()); 
				unsigned int j;
				j = mesh_.property(vparam_index_,vh_to); //Get Second Vertex in the Edge
				Mesh::EdgeHandle eh = mesh_.edge_handle(voh_it.handle());
				float weight = mesh_.property(eweight_,eh) / 2; //define weight as 1 in the edge
				W(i,j) = weight;
				total_weight += weight;
			}
		W(i,i) = -total_weight;
		}
	}


	//solve linear systems
	gmmMatrix W_t(W);
	solve_linear_system (W,by,y);
	solve_linear_system (W_t,bx,x);
}





void ParamViewer::
calc_harmonic_parameterization()
{
  // ------------- IMPLEMENT HERE ---------
	// TASK 5.2 harmonic map computation:
	// Search and order boundary vertices
  // Compute boundary parameters
  // Solve the linear system for internal vertex parameters using solve_linear_system()
  // Store the parameters in the vparam_h_ vertex property
	// ------------- IMPLEMENT HERE ---------
	gmmVector a,b; 
	calc_harmonic_parameterization_aux();
	calc_uniform_parameterization_aux(a,b);
	for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it)
	{
		int i = mesh_.property(vparam_index_,v_it.handle());
		Vec2f h(a[i],b[i]);
		mesh_.property(vparam_h_,v_it) = h;
    }

}

void ParamViewer::
calc_harmonic_parameterization_aux()
{
	
	Mesh::HalfedgeHandle h0, h1, h2;
	Mesh::VertexHandle vi, vj;
	Mesh::Point p0, p1, p2;
	Mesh::Point a,b;
	double w; 
    
	for (Mesh::EdgeIter e_it=mesh_.edges_begin(); e_it!=mesh_.edges_end(); ++e_it)
	{
		w  = 0.0;
		h0 = mesh_.halfedge_handle(e_it.handle(), 0);
		vi = mesh_.to_vertex_handle(h0);
		p0 = mesh_.point(vi);
		h1 = mesh_.halfedge_handle(e_it.handle(), 1);
		vj = mesh_.to_vertex_handle(h1);
		p1 = mesh_.point(vj);
		
		for (int i=0; i<2; i++)
		{
		Mesh::HalfedgeHandle h_choose;
			if (i == 0 ) h_choose = h0;
			else h_choose = h1;
			h2 = mesh_.next_halfedge_handle(h_choose);
			p2 = mesh_.point(mesh_.to_vertex_handle(h2));
			a = p0 - p2;
			b = p1 - p2;
			double cot_angle = dot (a,b) /  cross(a,b).length() ;
			w += cot_angle;
		}
		if (w < 0) w=0;

		mesh_.property(eweight_,e_it) = w;
	}
       
}
       

void 
ParamViewer::ComputeTextureCoordinates(int iTextureWidth, int iTextureHeight, int iRepeats, ParameterizationMode imode)
{
  // ------------- IMPLEMENT HERE ---------
	// TASK 5.3 Compute texture coordinates for textured mesh
  // rendering using a texture image of dimension iTextureWidth*iTextureHeights and iRepeats repeats
  // and store them in a mesh property.
  // If imode is equals to Uniform, compute the texture coordinates using the
  // parameters stored in vparam_u_ and store the result in texcoord_u_.
  // If imode is equals to Harmonic, compute the texture coordinates using the
  // parameters stored in vparam_h_ and store the result in texcoord_h_.
	// ------------- IMPLEMENT HERE ---------
	
	for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it)
	{
		Vec2f uv_param;
		if(imode == Uniform)
		{
			uv_param = mesh_.property(vparam_u_,v_it);
		}
		else if(imode == Harmonic)
		{
			uv_param = mesh_.property(vparam_h_,v_it);
		}	
		else
		{
			return;
		}
             
		Vec2f uv_n_repeats = (uv_param + Vec2f(1,1)) / 2;
		uv_n_repeats *= iRepeats;
		if(imode == Uniform)
		{
		 mesh_.property(texcoord_u_,v_it) = uv_n_repeats;
		}     
		else if(imode == Harmonic)
		{
			mesh_.property(texcoord_h_,v_it) = uv_n_repeats;
		}
		else
		{
			return;
		}
   }
}

//-----------------------------------------------------------------------------

void ParamViewer::calc_distortion(ParameterizationMode imode)
{
  float angle_distortion=0., area_distortion=0.;

  // ------------- IMPLEMENT HERE ---------
	// TASK 5.4 Compute distortion of triangle areas and angles
  // and print it in the output window.
  // If imode is equals to Uniform, uniform map distortion has to be 
  // computed.
  // If imode is equals to Harmonic, harmonic map distortion has to be
  // computed.
	// ------------- IMPLEMENT HERE ---------

  cout << "Parameterization distortion: " <<endl;
  cout << (imode==Uniform? " * Uniform map: ": " * Harmonic map: ") <<endl;
  cout << " ---- Angle distortion: " << angle_distortion << " Area distortion: " << area_distortion << endl;
}




bool
ParamViewer::
open_mesh(const char* _meshfilename)
{
	// load mesh
  if (MeshViewer::open_mesh(_meshfilename))
	{
    // store vertex initial positions and 3D mesh bounding box
    Mesh::VertexIter v_it=mesh_.vertices_begin(), v_end(mesh_.vertices_end());
    _bbMin3D = _bbMax3D = mesh_.point(v_it);
    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    {
      mesh_.property(vpos_,v_it) = mesh_.point(v_it);
      _bbMin3D.minimize(mesh_.point(v_it));
      _bbMax3D.maximize(mesh_.point(v_it));
	  _ParameterizationComputed_u = false;      
	  _ParameterizationComputed_h = false;      
	  if (_TextureCoordinates_u)
		delete [] _TextureCoordinates_u;

	  if (_TextureCoordinates_h)
		delete [] _TextureCoordinates_h;
	  _TextureCoordinates_u=NULL;
	  _TextureCoordinates_h=NULL;

    }
		return true;
	}
	return false;
}


void 
ParamViewer::
draw(const std::string& _draw_mode)
{
  if (indices_.empty())
  {
    MeshViewer::draw(_draw_mode);
    return;
  }
  if (_draw_mode == "UV Domain")
  {
    if (_ParameterizationMode_!=NoParameterization)
    {
      OpenMesh::VPropHandleT<OpenMesh::Vec2f> & TexCoordHandle = (_ParameterizationMode_ == Uniform? vparam_u_: vparam_h_);
      Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
      for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
      {
        OpenMesh::Vec2f UVCoord = mesh_.property(TexCoordHandle,v_it);
        mesh_.set_point(v_it, Mesh::Point(-UVCoord[0], UVCoord[1], 0.));
      }
      mesh_.update_normals();

      if (!_BoundingBox2DComputed)
      {
        Mesh::ConstVertexIter  v_it(mesh_.vertices_begin()), 
          v_end(mesh_.vertices_end());
        _bbMin2D = _bbMax2D = mesh_.point(v_it);
        for (; v_it!=v_end; ++v_it)
        {
          _bbMin2D.minimize(mesh_.point(v_it));
          _bbMax2D.maximize(mesh_.point(v_it));
        }
      }

      set_scene( (Vec3f)(_bbMin2D + _bbMax2D)*0.5, 0.5*(_bbMin2D - _bbMax2D).norm());

      MeshViewer::draw("Wireframe");
    }
  }
  else 
  {
    Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
      mesh_.set_point(v_it, mesh_.property(vpos_,v_it));
    mesh_.update_normals();

    set_scene( (Vec3f)(_bbMin3D + _bbMax3D)*0.5, 0.5*(_bbMin3D - _bbMax3D).norm());

    if (_draw_mode == "Textured mesh" && _ParameterizationMode_!=NoParameterization)
    {
      float *& TextureCoord = (_ParameterizationMode_ == Uniform? _TextureCoordinates_u: _TextureCoordinates_h);
      OpenMesh::VPropHandleT<OpenMesh::Vec2f> & TexCoordHandle = (_ParameterizationMode_ == Uniform? texcoord_u_: texcoord_h_);
      if (!TextureCoord)
      {
        int nvertices = mesh_.n_vertices();
        TextureCoord = new float[2*nvertices];
	  }


    Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
    int index=0;
    for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    {
      OpenMesh::Vec2f UVParam = mesh_.property(TexCoordHandle,v_it);
      TextureCoord[index++] = UVParam[0];
      TextureCoord[index++] = UVParam[1];
    }

      glEnable( GL_TEXTURE_2D ); 

      glEnable(GL_LIGHTING);
      glShadeModel(GL_SMOOTH);

      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_NORMAL_ARRAY);
      glEnableClientState(GL_TEXTURE_COORD_ARRAY);
      GL::glVertexPointer(mesh_.points());
      GL::glNormalPointer(mesh_.vertex_normals());
      GL::glTexCoordPointer(2, GL_FLOAT, 0, TextureCoord);

      glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

      glDisableClientState(GL_VERTEX_ARRAY);
      glDisableClientState(GL_NORMAL_ARRAY);
      glDisableClientState(GL_TEXTURE_COORD_ARRAY);

      glDisable( GL_TEXTURE_2D );
    }
    else MeshViewer::draw(_draw_mode);
  }
}


//-----------------------------------------------------------------------------

void
ParamViewer::
keyboard(int key, int x, int y)
{
  switch (toupper(key))
  { 
  case 'O':
	  {
		  //.pts (point cloud) file opening
		CFileDialog dlg(TRUE, LPCTSTR("off"), LPCTSTR("*.off"));
		if (dlg.DoModal() == IDOK){
			mesh_.clear();
			open_mesh(dlg.GetPathName().GetBuffer());
		}
	  break;
	  }
  case 'U':
    {
      _ParameterizationMode_ = Uniform;
      if (!_ParameterizationComputed_u)
      {
        calc_uniform_parameterization();
        ComputeTextureCoordinates(256, 256, _Repeat, _ParameterizationMode_);
        calc_distortion(_ParameterizationMode_);

        _ParameterizationComputed_u = true;      
      }
      glutPostRedisplay();
      break;
    }
  case 'H':
    {
      _ParameterizationMode_ = Harmonic;
      if (!_ParameterizationComputed_h)
      {
        calc_harmonic_parameterization();
        ComputeTextureCoordinates(256, 256, _Repeat, _ParameterizationMode_);
        calc_distortion(_ParameterizationMode_);

        _ParameterizationComputed_h = true;       
      }
      glutPostRedisplay();
      break;
    }
  case 'R': {
	  _Repeat++;
	  ComputeTextureCoordinates(256, 256, _Repeat, Uniform);
	  ComputeTextureCoordinates(256, 256, _Repeat, Harmonic);
	  printf("Number of repeats: %d\n",_Repeat);
			glutPostRedisplay();
			break;
			}
  case 'E': {
	  if (_Repeat>1) _Repeat--;
	  ComputeTextureCoordinates(256, 256, _Repeat, Uniform);
	  ComputeTextureCoordinates(256, 256, _Repeat, Harmonic);
	  printf("Number of repeats: %d\n",_Repeat);
			glutPostRedisplay();
			break;
		}

  default:
    {
      MeshViewer::keyboard(key, x, y);
      break;
    }
	}
}
//-----------------------------------------------------------------------------

//=============================================================================
