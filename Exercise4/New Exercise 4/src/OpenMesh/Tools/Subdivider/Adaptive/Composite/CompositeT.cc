//=============================================================================
//                                                                            
//                               OpenMesh                                     
//      Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen       
//                           www.openmesh.org                                 
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//                                License                                     
//                                                                            
//   This library is free software; you can redistribute it and/or modify it 
//   under the terms of the GNU Library General Public License as published  
//   by the Free Software Foundation, version 2.                             
//                                                                             
//   This library is distributed in the hope that it will be useful, but       
//   WITHOUT ANY WARRANTY; without even the implied warranty of                
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         
//   Library General Public License for more details.                          
//                                                                            
//   You should have received a copy of the GNU Library General Public         
//   License along with this library; if not, write to the Free Software       
//   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//   $Revision: 1.3 $
//   $Date: 2005/04/18 09:07:11 $
//                                                                            
//=============================================================================

/** \file Adaptive/Composite/CompositeT.cc
    
 */

//=============================================================================
//
//  CLASS CompositeT - IMPLEMENTATION
//
//=============================================================================

#define OPENMESH_SUBDIVIDER_ADAPTIVE_COMPOSITET_CC


//== INCLUDES =================================================================

#include <OpenMesh/Core/System/omstream.hh>
#include <OpenMesh/Core/System/config.hh>
#include <OpenMesh/Tools/Subdivider/Adaptive/Composite/CompositeT.hh>
#include <OpenMesh/Tools/Subdivider/Adaptive/Composite/RuleInterfaceT.hh>


//== NAMESPACE ================================================================

namespace OpenMesh   { // BEGIN_NS_OPENMESH
namespace Subdivider { // BEGIN_NS_DECIMATER
namespace Adaptive   { // BEGIN_NS_UNIFORM


//== IMPLEMENTATION ========================================================== 


template<class M>
bool
CompositeT<M> ::
initialize( void ) 
{
  typename Mesh::VertexIter  v_it;
  typename Mesh::FaceIter    f_it;
  typename Mesh::EdgeIter    e_it;
  const typename Mesh::Point zero_point(0.0, 0.0, 0.0);
  
  // ---------------------------------------- Init Vertices
  for (v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it) 
  {
    v_it->set_state(0);
    v_it->set_final();
    v_it->set_position(0, mesh_.point(v_it.handle()));
  }
  
  // ---------------------------------------- Init Faces
  for (f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); ++f_it) 
  {
    f_it->set_state(0);
    f_it->set_final();
    f_it->set_position(0, zero_point);
  }
  
  // ---------------------------------------- Init Edges
  for (e_it = mesh_.edges_begin(); e_it != mesh_.edges_end(); ++e_it) 
  {
    e_it->set_state(0);
    e_it->set_final();
    e_it->set_position(0, zero_point);
  }
  
  
  // ---------------------------------------- Init Rules

  int n_subdiv_rules_ = 0;


  // look for subdivision rule(s)
  for (size_t i=0; i < n_rules(); ++i) {

    if (rule_sequence_[i]->type()[0] == 'T' || 
        rule_sequence_[i]->type()[0] == 't') 
    {     
      ++n_subdiv_rules_;
      subdiv_rule_ = rule_sequence_[i];
      subdiv_type_ = rule_sequence_[i]->subdiv_type();
    }
  }


  // check for correct number of subdivision rules
  assert(n_subdiv_rules_ == 1);

  if (n_subdiv_rules_ != 1)
  {    
    std::cerr << "Error! More than one subdivision rules not allowed!\n";
    return false;
  }

  // check for subdivision type
  assert(subdiv_type_ == 3 || subdiv_type_ == 4);

  if (subdiv_type_ != 3 && subdiv_type_ != 4)
  {
    omerr() << "Error! Unknown subdivision type in sequence!" << std::endl;
    return false;
  }

  // set pointer to last rule
//   first_rule_ = rule_sequence_.front();
//   last_rule_  = rule_sequence_.back(); //[n_rules() - 1];

  // set numbers and previous rule
  for (size_t i = 0; i < n_rules(); ++i) 
  {
    rule_sequence_[i]->set_subdiv_type(subdiv_type_);
    rule_sequence_[i]->set_n_rules(n_rules());
    rule_sequence_[i]->set_number(i);
    rule_sequence_[i]->set_prev_rule(rule_sequence_[(i+n_rules()-1)%n_rules()]);
    rule_sequence_[i]->set_subdiv_rule(subdiv_rule_);
  }

  return true;
}


// ----------------------------------------------------------------------------
#define MOBJ mesh_.deref
#define TVH  to_vertex_handle
#define HEH  halfedge_handle
#define NHEH next_halfedge_handle
#define PHEH prev_halfedge_handle
#define OHEH opposite_halfedge_handle
// ----------------------------------------------------------------------------


template<class M>
void CompositeT<M>::refine(typename Mesh::FaceHandle& _fh) 
{
  std::vector<typename Mesh::HalfedgeHandle> hh_vector;

  // -------------------- calculate new level for faces and vertices
  int new_face_level = 
    t_rule()->number() + 1 + 
    ((int)floor((float)(MOBJ(_fh).state() - t_rule()->number() - 1)/n_rules()) + 1) * n_rules();

  int new_vertex_level = 
    new_face_level + l_rule()->number() - t_rule()->number();

  // -------------------- store old vertices
  // !!! only triangle meshes supported!
  typename Mesh::VertexHandle vh[3];

  vh[0] = mesh_.TVH(mesh_.HEH(_fh));
  vh[1] = mesh_.TVH(mesh_.NHEH(mesh_.HEH(_fh)));
  vh[2] = mesh_.TVH(mesh_.PHEH(mesh_.HEH(_fh)));

  // save handles to incoming halfedges for getting the new vertices 
  // after subdivision (1-4 split)
  if (subdiv_type_ == 4) 
  {
    hh_vector.clear();

    // green face
    if (MOBJ(_fh).final()) 
    {
      typename Mesh::FaceHalfedgeIter fh_it(mesh_.fh_iter(_fh));

      for (; fh_it; ++fh_it) 
      {
	hh_vector.push_back(mesh_.PHEH(mesh_.OHEH(fh_it.handle())));
      }
    }

    // red face
    else 
    {
      
      typename Mesh::HalfedgeHandle red_hh(MOBJ(_fh).red_halfedge());

      hh_vector.push_back(mesh_.PHEH(mesh_.OHEH(mesh_.NHEH(red_hh))));
      hh_vector.push_back(mesh_.PHEH(mesh_.OHEH(mesh_.PHEH(mesh_.OHEH(red_hh)))));
    }
  }

  
  // -------------------- Average rule before topo rule?
  if (t_rule()->number() > 0)
    t_rule()->prev_rule()->raise(_fh, new_face_level-1);

  // -------------------- Apply topological operator first
  t_rule()->raise(_fh, new_face_level);

#if 0 // original code
  assert(MOBJ(_fh).state() >= 
         subdiv_rule_->number()+1+(int) (MOBJ(_fh).state()/n_rules())*n_rules());
#else // improved code (use % operation and avoid floating point division)
  assert( MOBJ(_fh).state() >= ( t_rule()->number()+1+generation(_fh) ) );
#endif

  // raise new vertices to final levels
  if (subdiv_type_ == 3) 
  {
    typename Mesh::VertexHandle new_vh(mesh_.TVH(mesh_.NHEH(mesh_.HEH(_fh))));

    // raise new vertex to final level
    l_rule()->raise(new_vh, new_vertex_level);
  }

  if (subdiv_type_ == 4) 
  {
    typename Mesh::HalfedgeHandle hh;
    typename Mesh::VertexHandle   new_vh;

    while (!hh_vector.empty()) {

      hh = hh_vector.back();
      hh_vector.pop_back();

      // get new vertex
      new_vh = mesh_.TVH(mesh_.NHEH(hh));

      // raise new vertex to final level
      l_rule()->raise(new_vh, new_vertex_level);
    }
  }

  // raise old vertices to final position 
  l_rule()->raise(vh[0], new_vertex_level);
  l_rule()->raise(vh[1], new_vertex_level);
  l_rule()->raise(vh[2], new_vertex_level);
}


// ----------------------------------------------------------------------------


template<class M>
void CompositeT<M>::refine(typename Mesh::VertexHandle& _vh) 
{
  // calculate next final level for vertex
  int new_vertex_state = generation(_vh) + l_rule()->number() + 1;

  assert( new_vertex_state == mesh_.deref(_vh).state()+1 );

  // raise vertex to final position
  l_rule()->raise(_vh, new_vertex_state);
}


// ----------------------------------------------------------------------------


template <class M>
std::string CompositeT<M>::rules_as_string(const std::string& _sep) const
{
  std::string seq;
  typename RuleSequence::const_iterator it = rule_sequence_.begin();

  if ( it != rule_sequence_.end() )
  {
    seq = (*it)->type();
    for (++it; it != rule_sequence_.end(); ++it )
    {
      seq += _sep;
      seq += (*it)->type();    
    }
  }
  return seq;
}

// ----------------------------------------------------------------------------
#undef MOBJ
#undef TVH
#undef HEH
#undef NHEH
#undef PHEH
#undef OHEH
//=============================================================================
} // END_NS_ADAPTIVE
} // END_NS_SUBDIVIDER
} // END_NS_OPENMESH
//=============================================================================
