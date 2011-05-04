//=============================================================================
//                                                
//   Code framework for the lecture
//
//   "Surface Representation and Geometric Modeling"
//
//   Mark Pauly, Mario Botsch, Balint Miklos, and Hao Li
//
//   Copyright (C) 2007 by  Computer Graphics Laboratory, ETH Zurich
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
//  CLASS MeshViewerWidget
//
//=============================================================================


#ifndef MESH_VIEWER_WIDGET_HH
#define MESH_VIEWER_WIDGET_HH


//== INCLUDES =================================================================


#include "GlutExaminer.hh"
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include "QuadricT.hh"
#include <set>
#include <fstream>





//== CLASS DEFINITION =========================================================


class MeshViewer : public GlutExaminer
{

protected:
	//std::fstream out;
  typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;
  int percentage_;
  
public:
   
  /// default constructor
  MeshViewer(const char* _title, int _width, int _height);

  /// open mesh
  virtual bool open_mesh(const char* _filename);

  /// update buffer with face indices
  void update_face_indices();

  /// draw the scene
  virtual void draw(const std::string& _draw_mode);

  virtual void MeshViewer::keyboard(int key, int x, int y);

  //mesh decimation actions
  void  init();

  void CalculateVertexQuadric( Mesh::VertexHandle vh );
  void CalculateVertexPriority( Mesh::VertexHandle vh );

  bool  is_collapse_legal(Mesh::HalfedgeHandle _hh);
  float priority(Mesh::HalfedgeHandle _heh);
  void  decimate(unsigned int _n_vertices);
  void  enqueue_vertex(Mesh::VertexHandle vh);
  void  dequeue_vertex(Mesh::VertexHandle vh);

  // access quadric of vertex _vh
	Quadricd& quadric(Mesh::VertexHandle _vh) 
	{ return mesh_.property(vquadric, _vh); }

	// access priority of vertex _vh
	float& priority(Mesh::VertexHandle _vh)
	{ return mesh_.property(vprio, _vh); }

	// access target halfedge of vertex _vh
	Mesh::HalfedgeHandle& target(Mesh::VertexHandle _vh)
	{ return mesh_.property(vtarget, _vh); }

protected:

	  Mesh                       mesh_;
  std::vector<unsigned int>  indices_;

  OpenMesh::VPropHandleT<Quadricd>                vquadric;
  OpenMesh::VPropHandleT<float>                   vprio;
  OpenMesh::VPropHandleT<Mesh::HalfedgeHandle>  vtarget;
  OpenMesh::HPropHandleT<bool>		vmydeletedstatus;

  struct QueueVertex{
	  Mesh::VertexHandle v;
	  float prio;
  };

  // compare functor for priority queue
  struct VertexCmp
	{

		bool operator()(QueueVertex _v0, QueueVertex _v1) const
		{
			// std::set needs UNIQUE keys -> handle equal priorities
			return (( _v0.prio ==  _v1.prio) ? 
				(_v0.v.idx() < _v1.v.idx()) :
			(  _v0.prio <   _v1.prio));
		}
	};


	std::set<QueueVertex, VertexCmp>  queue;

	Mesh::Normal CalculateNormal( VertexHandle a, VertexHandle b, VertexHandle c)
	{	
		Mesh::Point pa, pb, pc;
		pa = mesh_.point(a);
		pb = mesh_.point(b);
		pc = mesh_.point(c);
		return CalculateNormal(pa,pb,pc);
	}

	Mesh::Normal CalculateNormal( VertexHandle a, VertexHandle b, VertexHandle c, Mesh::VertexHandle v0, Mesh::Point p1 )
	{
		Mesh::Point pa, pb, pc;
		pa = mesh_.point(a);
		pb = mesh_.point(b);
		pc = mesh_.point(c);
		if (a == v0) return CalculateNormal(p1,pb,pc);
		else if (b == v0) return CalculateNormal(pa,p1,pc);
		else return CalculateNormal(pa,pb,p1);	
	}

	Mesh::Normal CalculateNormal( Mesh::Point a, Mesh::Point b, Mesh::Point c)
	{
		OpenMesh::Vec3f v1 = b - a;
		OpenMesh::Vec3f v2 = c - a;
		return cross(v1,v2).normalize();
	}
};


//=============================================================================
#endif // MESH_VIEWER_WIDGET_HH defined
//=============================================================================

