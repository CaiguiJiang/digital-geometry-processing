//=============================================================================
//                                                                            
//   Example code for the full-day course
//
//   M. Botsch, M. Pauly, C. Roessl, S. Bischoff, L. Kobbelt,
//   "Geometric Modeling Based on Triangle Meshes"
//   held at SIGGRAPH 2006, Boston, and Eurographics 2006, Vienna.
//
//   Copyright (C) 2006 by  Computer Graphics Laboratory, ETH Zurich, 
//                      and Computer Graphics Group,      RWTH Aachen
//
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//                                License                                     
//                                                                            
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//   
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//   
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor, 
//   Boston, MA  02110-1301, USA.
//                                                                            
//=============================================================================
//=============================================================================
//
//  CLASS SmoothingViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "stdafx.h"
#include "SmoothingViewer.hh"
#include <vector>


//== IMPLEMENTATION ========================================================== 


SmoothingViewer::
SmoothingViewer(const char* _title, int _width, int _height)
  : QualityViewer(_title, _width, _height)
{ 
  mesh_.add_property(vpos_);
}


//-----------------------------------------------------------------------------


void
SmoothingViewer::
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
	case 'N':
		{
			std::cout << "10 Laplace-Beltrami smoothing iterations: " << std::flush;
			smooth(10);
			calc_weights();
			calc_mean_curvature();
			calc_uniform_mean_curvature();
			calc_gauss_curvature();
			calc_triangle_quality();
			face_color_coding();

			glutPostRedisplay();
			std::cout << "done\n";
			break;
		}
	case 'U':
		{
			std::cout << "10 uniform smoothing iterations: " << std::flush;
			uniform_smooth(10);
			calc_weights();
			calc_mean_curvature();
			calc_uniform_mean_curvature();
			calc_gauss_curvature();
			calc_triangle_quality();
			face_color_coding();

			glutPostRedisplay();
			std::cout << "done\n";
			break;
		}


	default:
		{
			QualityViewer::keyboard(key, x, y);
			break;
		}
	}
}


//-----------------------------------------------------------------------------


void 
SmoothingViewer::
smooth(unsigned int _iters)
{
	// ------------- IMPLEMENT HERE ---------
	// TASK 3.3.b Smoothing using the Laplace-Beltrami.
	// Use eweight_ properties for the individual edge weights
	// and their sum for the normalization term.
	// ------------- IMPLEMENT HERE ---------



}

//-----------------------------------------------------------------------------


void 
SmoothingViewer::
uniform_smooth(unsigned int _iters)
{
	// ------------- IMPLEMENT HERE ---------
	// TASK 3.1.b Smoothing using the uniform Laplacian approximation
	// ------------- IMPLEMENT HERE ---------

	

	for (unsigned int i = 0; i < _iters; i++) {
		std::vector<Mesh::Point> Lu;
		for (Mesh::VertexIter vIt = mesh_.vertices_begin();
		vIt != mesh_.vertices_end(); ++vIt)
		{
			int valence = 0;
			Mesh::Point Lu_v;
			Lu_v[0] = Lu_v[1] = Lu_v[2] = 0.0;
			
			// Calculate number of neighbors around current vertex, and vector sum of all neighbors
			for (Mesh::VertexVertexIter vvIt = mesh_.vv_iter(vIt.handle());	vvIt; ++vvIt, ++valence) {
				Lu_v += mesh_.point(vvIt.handle());
			}

			Lu_v /= valence;
			Lu_v -= mesh_.point(vIt.handle());
			Lu.push_back(Lu_v);
		}
		
		int i = 0;
		for (Mesh::VertexIter vIt = mesh_.vertices_begin();
			 vIt != mesh_.vertices_end(); ++vIt, ++i) {
				 Mesh::Point newV = mesh_.point(vIt.handle()) + Lu[i] * 0.5;
				 mesh_.set_point(vIt.handle(), newV);
		}
		mesh_.update_normals();
	}
}

//=============================================================================
