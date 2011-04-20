//=============================================================================
//                                                
//   Code framework for the lecture
//
//   "Surface Representation and Geometric Modeling"
//
//   Mark Pauly, Mario Botsch, Balint Miklos, and Hao Li
//
//   Copyright (C) 2007 by  Applied Geometry Group and 
//							Computer Graphics Laboratory, ETH Zurich
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
//  CLASS MeshViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================


#include "StdAfx.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshViewer.hh"
#include "gl.hh"
#include <iostream>
#include <fstream>
#include <set>
//#include <unistd.h> //in linux version this header file was included. Removed for MSVC version
#include <float.h>
#include "OpenMesh/Tools/Geometry/QuadricT.hh"


//== IMPLEMENTATION ========================================================== 


MeshViewer::
MeshViewer(const char* _title, int _width, int _height)
  : GlutExaminer(_title, _width, _height)
{
  mesh_.request_face_normals();
  mesh_.request_vertex_normals();

  clear_draw_modes();
  add_draw_mode("Wireframe");
  add_draw_mode("Hidden Line");
  add_draw_mode("Solid Flat");
  add_draw_mode("Solid Smooth");
  set_draw_mode(3);

  // add required properties for decimation
  mesh_.request_vertex_status();
  mesh_.request_edge_status();
  mesh_.request_face_status();
  mesh_.add_property(vquadric);
  mesh_.add_property(vprio);
  mesh_.add_property(vtarget);
  percentage_=50;
}


//-----------------------------------------------------------------------------


bool
MeshViewer::
open_mesh(const char* _filename)
{
  // load mesh
  if (OpenMesh::IO::read_mesh(mesh_, _filename))
  {
    // set center and radius
    Mesh::ConstVertexIter  v_it(mesh_.vertices_begin()), 
                           v_end(mesh_.vertices_end());
    Mesh::Point            bbMin, bbMax;

    bbMin = bbMax = mesh_.point(v_it);
    for (; v_it!=v_end; ++v_it)
    {
      bbMin.minimize(mesh_.point(v_it));
      bbMax.maximize(mesh_.point(v_it));
    }
    set_scene( (Vec3f)(bbMin + bbMax)*0.5, 0.5*(bbMin - bbMax).norm());


    // compute face & vertex normals
    mesh_.update_normals();


    // update face indices for faster rendering
    update_face_indices();

    // info
    std::cerr << mesh_.n_vertices() << " vertices, "
	      << mesh_.n_faces()    << " faces\n";

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------------


void
MeshViewer::
update_face_indices()
{
  Mesh::ConstFaceIter        f_it(mesh_.faces_sbegin()), 
                             f_end(mesh_.faces_end());
  Mesh::ConstFaceVertexIter  fv_it;

  indices_.clear();
  indices_.reserve(mesh_.n_faces()*3);
  std::cout << "mesh indices updated" << std::endl;

  for (; f_it!=f_end; ++f_it)
  {
    indices_.push_back((fv_it=mesh_.cfv_iter(f_it)).handle().idx());
    indices_.push_back((++fv_it).handle().idx());
    indices_.push_back((++fv_it).handle().idx());
  }
}


//-----------------------------------------------------------------------------


void 
MeshViewer::
draw(const std::string& _draw_mode)
{
  if (indices_.empty())
  {
    GlutExaminer::draw(_draw_mode);
    return;
  }



  if (_draw_mode == "Wireframe")
  {
    glDisable(GL_LIGHTING);
    glColor3f(1.0, 1.0, 1.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glEnableClientState(GL_VERTEX_ARRAY);
    GL::glVertexPointer(mesh_.points());

    glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

    glDisableClientState(GL_VERTEX_ARRAY);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }


  else if (_draw_mode == "Hidden Line")
  {

	  glDisable(GL_LIGHTING);
	  glShadeModel(GL_SMOOTH);
	  glColor3f(0.0, 0.0, 0.0);

	  glEnableClientState(GL_VERTEX_ARRAY);
	  GL::glVertexPointer(mesh_.points());

	  glDepthRange(0.01, 1.0);
	  glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);
	  glDisableClientState(GL_VERTEX_ARRAY);
	  glColor3f(1.0, 1.0, 1.0);

	  glEnableClientState(GL_VERTEX_ARRAY);
	  GL::glVertexPointer(mesh_.points());

	  glDrawBuffer(GL_BACK);
	  glDepthRange(0.0, 1.0);
	  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	  glDepthFunc(GL_LEQUAL);
	  glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

	  glDisableClientState(GL_VERTEX_ARRAY);
	  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	  glDepthFunc(GL_LESS);

  }


  else if (_draw_mode == "Solid Flat")
  {
    Mesh::ConstFaceIter        f_it(mesh_.faces_begin()), 
                               f_end(mesh_.faces_end());
    Mesh::ConstFaceVertexIter  fv_it;

    glEnable(GL_LIGHTING);
    glShadeModel(GL_FLAT);

    glBegin(GL_TRIANGLES);
    for (; f_it!=f_end; ++f_it)
    {
      GL::glNormal(mesh_.normal(f_it));
      fv_it = mesh_.cfv_iter(f_it.handle()); 
      GL::glVertex(mesh_.point(fv_it));
      ++fv_it;
      GL::glVertex(mesh_.point(fv_it));
      ++fv_it;
      GL::glVertex(mesh_.point(fv_it));
    }
    glEnd();
  }


  else if (_draw_mode == "Solid Smooth")
  {
    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    GL::glVertexPointer(mesh_.points());
    GL::glNormalPointer(mesh_.vertex_normals());

    glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
  }
}

void
MeshViewer::
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
		case 'Z':  //up percentage by 5%
			{
				percentage_=percentage_+5;
				if (percentage_>95) percentage_=95;
				std::cout<<"Percentage is %"<<percentage_<<"\n";
			}
			break;

		case 'X':  //up percentage by 5%
			{
				percentage_=percentage_-5;
				if (percentage_<5) percentage_=5;
				std::cout<<"Percentage is %"<<percentage_<<"\n";
			}
			break;

		case 'D': //decimate
			{
				// compute normals & quadrics
				init();

				// decimate
				decimate(percentage_ * mesh_.n_vertices() / 100);
				std::cout << "#vertices: " << mesh_.n_vertices() << std::endl;
				break;
			}
		case 'M':
			{
				// compute normals & quadrics
				init();

				// decimate
				decimate(mesh_.n_vertices() - 1);
				std::cout << "#vertices: " << mesh_.n_vertices() << std::endl;
				break;
			}

		default:
			{
				GlutExaminer::keyboard(key, x, y);
				break;
			}
	}
}


//=============================================================================
//DECIMATION IMPLEMENTATION FUNCTIONS
//==============================================================================




void MeshViewer::init()
{
	// compute face normals
	mesh_.update_face_normals();

	Mesh::VertexIter  v_it, v_end = mesh_.vertices_end();	

	for (v_it=mesh_.vertices_begin(); v_it != v_end; ++v_it)
	{
		priority(v_it) = -1.0;
		quadric(v_it).clear();
		
		Mesh::Point v = mesh_.point(v_it.handle());

		// iterate through all adjacent faces to current vertex
		for (Mesh::VertexFaceIter vf_it = mesh_.vf_iter(v_it.handle()); vf_it; ++vf_it)
		{	
			Mesh::Normal n = mesh_.normal(vf_it.handle());

			double a = n[0];
			double b = n[1];
			double c = n[2];
			double d = -n[0]*v[0]-n[1]*v[1]-n[2]*v[2];						
			Quadricd q(a,b,c,d);
			
			quadric(v_it) += q;			
		}
	}
}


//-----------------------------------------------------------------------------


bool MeshViewer::is_collapse_legal(Mesh::HalfedgeHandle _hh)
{
	// collect vertices
	Mesh::VertexHandle v0, v1;
	v0 = mesh_.from_vertex_handle(_hh);
	v1 = mesh_.to_vertex_handle(_hh);

	// collect faces
	Mesh::FaceHandle fl = mesh_.face_handle(_hh);
	Mesh::FaceHandle fr = mesh_.face_handle(mesh_.opposite_halfedge_handle(_hh));

	// backup point positions
	Mesh::Point p0 = mesh_.point(v0);
	Mesh::Point p1 = mesh_.point(v1);

	// topological test
	if (!mesh_.is_collapse_ok(_hh))
		return false;

	// test boundary stuff
	if (mesh_.is_boundary(v0) && !mesh_.is_boundary(v1))
		return false;

	// Exercise 4.2 -----------------------------------------------
	// INSERT CODE:
	// test normal flipping:
	//   if normal vector of a (non-degenerate) triangle changes by 
	//   more than pi/4 degrees, return false.
	// ------------------------------------------------------------

	// Iterate through all faces adjacent to source vertex
	bool legal;
	for (Mesh::VertexFaceIter vf_it = mesh_.vf_iter(v0); vf_it; ++vf_it)
	{		
		Mesh::FaceHandle fh = vf_it.handle();
		// Skip soon-to-be degenerate triangles
		if (fh == fl || fh == fr) continue;
		
		// Iterate through		
		Mesh::ConstFaceVertexIter cfv_it = mesh_.cfv_iter(fh);
		VertexHandle a = cfv_it.handle(); ++cfv_it;
		VertexHandle b = cfv_it.handle(); ++cfv_it;
		VertexHandle c = cfv_it.handle();
		
		Mesh::Normal n_old = CalculateNormal(a,b,c,v0);
		Mesh::Normal n_new = CalculateNormal(a,b,c,v0,p1);

		double dot_product = dot(n_old, n_new);
		if (acos(dot_product) > 0.78539816339744830961566084581988) 
		{
			legal = false;
			break;
		}
	}

	// collapse passed all tests -> ok
	return legal;
}


//-----------------------------------------------------------------------------


float MeshViewer::priority(Mesh::HalfedgeHandle _heh)
{
	Mesh::VertexHandle vh0, vh1;
	vh0 = mesh_.from_vertex_handle(_heh);
	vh1 = mesh_.to_vertex_handle(_heh);
		
	Quadricd& d0 = quadric(vh0);
	Quadricd& d1 = quadric(vh1);
	
	OpenMesh::Vec3d v0,v1;
	v0 = mesh_.point(vh0);
	v1 = mesh_.point(vh1);

	return d0(v0) + d1(v1);
}


//-----------------------------------------------------------------------------


void MeshViewer::enqueue_vertex(Mesh::VertexHandle _vh)
{
	float                   prio, min_prio(FLT_MAX);
	Mesh::HalfedgeHandle  min_hh;


	// find best out-going halfedge
	for (Mesh::VOHIter vh_it(mesh_, _vh); vh_it; ++vh_it)
	{
		if (is_collapse_legal(vh_it))
		{
			prio = priority(vh_it);
			if (prio >= -1.0 && prio < min_prio)
			{
				min_prio = prio;
				min_hh   = vh_it.handle();
			}
		}
	}


	// update queue
	QueueVertex qv;
	qv.v=_vh; qv.prio=priority(_vh);
	if (priority(_vh) != -1.0) 
	{
		queue.erase(qv);
		priority(_vh) = -1.0;
	}

	if (min_hh.is_valid()) 
	{
		priority(_vh) = min_prio;
		target(_vh)   = min_hh;
		qv.prio=min_prio;
		queue.insert(qv);
	}
}


//-----------------------------------------------------------------------------


void MeshViewer::decimate(unsigned int _n_vertices)
{
	unsigned int nv(mesh_.n_vertices());

	Mesh::HalfedgeHandle hh;
	Mesh::VertexHandle   to, from;
	Mesh::VVIter         vv_it;

	std::vector<Mesh::VertexHandle>            one_ring;
	std::vector<Mesh::VertexHandle>::iterator  or_it, or_end;



	// build priority queue
	Mesh::VertexIter  v_it  = mesh_.vertices_begin(), 
		v_end = mesh_.vertices_end();

	queue.clear();
	for (; v_it!=v_end; ++v_it)
		enqueue_vertex(v_it.handle());

	while (nv > _n_vertices && !queue.empty())
	{
		//std::cout << "# Vertices reduced to " << nv << ". " << queue.size() << " vertices remain in queue." << std::endl;
		// Exercise 4.3 ----------------------------------------------
		// INSERT CODE:
		// Decimate using priority queue:
		//   1) take 1st element of queue
		//   2) collapse this halfedge
		//   3) update queue
		// -----------------------------------------------------------
		std::set<QueueVertex, VertexCmp>::iterator it = queue.begin();
		from = it->v;
		queue.erase(it);
		
		hh = target(from);
		//if (is_collapse_legal(hh))
		//{
			for (vv_it = mesh_.vv_iter(from); vv_it; ++vv_it)
			{
				one_ring.push_back(vv_it.handle());
			}


			mesh_.collapse(hh);
			
			for (or_it = one_ring.begin(); or_it != one_ring.end(); ++or_it)
			{
				enqueue_vertex(*or_it);
			}
		//}

		nv--;
	}



	// clean up
	queue.clear();

	// now, delete the items marked to be deleted
	mesh_.garbage_collection();

	// re-compute face & vertex normals
	mesh_.update_normals();


	// re-update face indices for faster rendering
	update_face_indices();
}
