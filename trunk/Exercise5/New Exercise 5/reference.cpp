
   8.
      //== INCLUDES =================================================================
   9.
       
  10.
      #include "HarmonicMapViewer.hh"
  11.
       
  12.
      //== IMPLEMENTATION ==========================================================
  13.
       
  14.
       
  15.
      HarmonicMapViewer::
  16.
      HarmonicMapViewer(const char* _title, int _width, int _height, int iRepeat)
  17.
        : MeshViewer(_title, _width, _height)
  18.
      {
  19.
        _Repeat = iRepeat;
  20.
       
  21.
        mesh_.request_vertex_colors();
  22.
       
  23.
        mesh_.add_property(vpos_);
  24.
        mesh_.add_property(vparam_u_);
  25.
        mesh_.add_property(vparam_h_);
  26.
        mesh_.add_property(texcoord_u_);
  27.
        mesh_.add_property(texcoord_h_);
  28.
        mesh_.add_property(vparam_index_);
  29.
        mesh_.add_property(eweight_);
  30.
       
  31.
        add_draw_mode("UV Domain");
  32.
              add_draw_mode("Textured mesh");
  33.
       
  34.
        _TextureCoordinates_u=NULL;
  35.
        _TextureCoordinates_h=NULL;
  36.
       
  37.
        _ParameterizationMode_ = NoParameterization;
  38.
       
  39.
              init();
  40.
      }
  41.
       
  42.
      HarmonicMapViewer::
  43.
      ~HarmonicMapViewer()
  44.
      {
  45.
        if (glIsTexture(textureID_))  
  46.
                      glDeleteTextures( 1, &textureID_);
  47.
       
  48.
        if (_TextureCoordinates_u)
  49.
          delete [] _TextureCoordinates_u;
  50.
       
  51.
        if (_TextureCoordinates_h)
  52.
          delete [] _TextureCoordinates_h;
  53.
      }
  54.
       
  55.
      void
  56.
      HarmonicMapViewer::
  57.
      init()
  58.
      {
  59.
              // base class first
  60.
              MeshViewer::init();
  61.
       
  62.
        // generate checkerboard-like image
  63.
        GLubyte tex[256][256][3];
  64.
        int index=0;
  65.
        for (int x=0; x<256; ++x)
  66.
        {
  67.
          for (int y=0; y<256; ++y)
  68.
          {
  69.
            if ((x<128&&y<128) ||(x>128&&y>128))
  70.
            {
  71.
              tex[x][y][0] = 0;
  72.
              tex[x][y][1] = 255;
  73.
              tex[x][y][2] = 0;
  74.
            }
  75.
            else
  76.
            {
  77.
              tex[x][y][0] = 255;
  78.
              tex[x][y][1] = 255;
  79.
              tex[x][y][2] = 255;
  80.
            }
  81.
          }
  82.
        }
  83.
        // generate texture
  84.
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  85.
        if (!glIsTexture(textureID_))
  86.
          glGenTextures(1, &textureID_);
  87.
        glBindTexture(GL_TEXTURE_2D, textureID_);
  88.
       
  89.
        // copy texture to GL
  90.
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  91.
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  92.
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  93.
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  94.
        glTexImage2D(GL_TEXTURE_2D, 0, 3, 256, 256,
  95.
          0, GL_RGB, GL_UNSIGNED_BYTE, tex);
  96.
      }
  97.
       
  98.
      bool
  99.
      HarmonicMapViewer::
 100.
      open_mesh(const char* _meshfilename)
 101.
      {
 102.
              // load mesh
 103.
        if (MeshViewer::open_mesh(_meshfilename))
 104.
              {
 105.
          // store vertex initial positions and 3D mesh bounding box
 106.
          Mesh::VertexIter v_it=mesh_.vertices_begin(), v_end(mesh_.vertices_end());
 107.
          _bbMin3D = _bbMax3D = mesh_.point(v_it);
 108.
          for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
 109.
          {
 110.
            mesh_.property(vpos_,v_it) = mesh_.point(v_it);
 111.
            _bbMin3D.minimize(mesh_.point(v_it));
 112.
            _bbMax3D.maximize(mesh_.point(v_it));
 113.
          }
 114.
                      return true;
 115.
              }
 116.
              return false;
 117.
      }
 118.
       
 119.
      //-----------------------------------------------------------------------------
 120.
       
 121.
       
 122.
      void
 123.
      HarmonicMapViewer::solve_linear_system( gmmMatrix& _M, gmmVector& _b, gmmVector& _x)
 124.
      {
 125.
        unsigned int N = _b.size();
 126.
              _x.resize(N);
 127.
        std::vector< size_t >  ipvt(N);
 128.
        gmm::lu_factor( _M, ipvt );
 129.
        gmm::lu_solve( _M, ipvt, _x, _b );
 130.
      }
 131.
       
 132.
       
 133.
       
 134.
      void HarmonicMapViewer::
 135.
      calc_parameterization(gmmVector& x, gmmVector& y){
 136.
       
 137.
        // Init data
 138.
        size_t n = mesh_.n_vertices();
 139.
        gmmMatrix W(n,n);
 140.
        gmmVector bx(n),by(n);
 141.
        Mesh::Scalar PI = atan(1.)*4;
 142.
       
 143.
        // Find a border half edge
 144.
        HalfedgeHandle h_init_boundary;
 145.
        bool found_boundary = false;
 146.
        for(Mesh::HalfedgeIter h_it=mesh_.halfedges_begin(); h_it!=mesh_.halfedges_end(); ++h_it){
 147.
          if(mesh_.is_boundary(h_it)){
 148.
            h_init_boundary = h_it.handle();
 149.
            found_boundary = true;
 150.
            break;
 151.
          }
 152.
        }
 153.
        assert(found_boundary);
 154.
       
 155.
        // Calculate total boundary length
 156.
        Mesh::HalfedgeHandle h_cur = h_init_boundary;
 157.
        Mesh::Scalar boundary_length = 0;
 158.
        unsigned int n_boundary = 0;
 159.
        do {
 160.
          boundary_length += mesh_.calc_edge_length(h_cur);
 161.
          h_cur = mesh_.next_halfedge_handle(h_cur);
 162.
          n_boundary++;
 163.
        } while(h_cur != h_init_boundary);
 164.
       
 165.
       
 166.
        // Iterate over boundary
 167.
        h_cur = h_init_boundary;
 168.
        unsigned int i=0;
 169.
        Mesh::Scalar theta=0;
 170.
        do {
 171.
       
 172.
          // Get vertex, save vertex index in matrix
 173.
          Mesh::VertexHandle vh = mesh_.to_vertex_handle(h_cur);
 174.
          mesh_.property(vparam_index_,vh) = i;
 175.
       
 176.
          W(i,i) = 1;
 177.
          bx[i] = cos(theta);
 178.
          by[i] = sin(theta);
 179.
       
 180.
          // Next boundary halfedge
 181.
          h_cur = mesh_.next_halfedge_handle(h_cur);
 182.
          i++;
 183.
          theta += (mesh_.calc_edge_length(h_cur)/boundary_length)*2*PI;
 184.
        } while(h_cur != h_init_boundary);
 185.
       
 186.
        // Set the indices for all non-boundary vertices
 187.
        for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){
 188.
          if(!mesh_.is_boundary(v_it)){
 189.
            mesh_.property(vparam_index_,v_it.handle()) = i++;
 190.
          }      
 191.
        }
 192.
       
 193.
        // Set matrix values for non-boundary vertices
 194.
        for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){
 195.
          if(!mesh_.is_boundary(v_it)){
 196.
            unsigned int i = mesh_.property(vparam_index_,v_it);
 197.
            Mesh::Scalar total_weight = 0;
 198.
       
 199.
            // All outgoing half edges
 200.
            for(Mesh::VertexOHalfedgeIter voh_it=mesh_.voh_iter(v_it.handle()); voh_it;++voh_it){
 201.
              Mesh::VertexHandle vh_to = mesh_.to_vertex_handle(voh_it.handle());
 202.
              unsigned int j = mesh_.property(vparam_index_,vh_to);
 203.
              Mesh::EdgeHandle eh = mesh_.edge_handle(voh_it.handle());
 204.
              Mesh::Scalar weight = mesh_.property(eweight_,eh) / 2; // XXX: Notice division!
 205.
              W(i,j) = weight;
 206.
              total_weight += weight;
 207.
            }
 208.
       
 209.
            W(i,i) = -total_weight;
 210.
          }      
 211.
        }
 212.
       
 213.
        //for(unsigned int i=0;i<n;i++){
 214.
        //  for(unsigned int j=0;j<n;j++)
 215.
        //    cout << W(i,j) << " ";
 216.
        //  
 217.
        //  cout << endl;
 218.
        //}
 219.
       
 220.
        // Solve linear systems
 221.
        gmmMatrix W_t(W);
 222.
        solve_linear_system(W,by,y);
 223.
        solve_linear_system(W_t,bx,x);
 224.
       
 225.
      }
 226.
       
 227.
      //-----------------------------------------------------------------------------
 228.
       
 229.
      void HarmonicMapViewer::
 230.
      calc_uniform_parameterization()
 231.
      {
 232.
        // ------------- IMPLEMENT HERE ---------
 233.
              // TASK 5.1 Uniform map computation:
 234.
              // Search and order boundary vertices
 235.
        // Compute boundary parameters
 236.
        // Solve the linear system for internal vertex parameters using solve_linear_system()
 237.
        // Store the parameters in the vparam_u_ vertex property
 238.
              // ------------- IMPLEMENT HERE ---------
 239.
       
 240.
        // Init weights
 241.
        for(Mesh::EdgeIter e_it=mesh_.edges_begin();e_it!=mesh_.edges_end();++e_it){
 242.
          // Weight will be divided evenly between vertices, hence resulting vertice weights will be 1!
 243.
          mesh_.property(eweight_, e_it) = 2;
 244.
        }
 245.
       
 246.
        // Calculate parameterization based on eweights_ edge properties
 247.
        gmmVector x,y;
 248.
        calc_parameterization(x,y);
 249.
       
 250.
        // Save to vertices
 251.
        for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){
 252.
          unsigned int i = mesh_.property(vparam_index_,v_it.handle());
 253.
          Vec2f u(x[i],y[i]);
 254.
          mesh_.property(vparam_u_,v_it) = u;
 255.
        }
 256.
       
 257.
      }
 258.
       
 259.
      void
 260.
      HarmonicMapViewer::
 261.
      calc_weights()
 262.
      {
 263.
              Mesh::VertexIter        v_it, v_end(mesh_.vertices_end());
 264.
              Mesh::EdgeIter          e_it, e_end(mesh_.edges_end());
 265.
              Mesh::VertexFaceIter    vf_it;
 266.
              Mesh::FaceVertexIter    fv_it;
 267.
              Mesh::HalfedgeHandle    h0, h1, h2;
 268.
              Mesh::VertexHandle      v0, v1;
 269.
              Mesh::Point             p0, p1, p2, d0, d1;
 270.
              Mesh::Scalar            w;
 271.
       
 272.
       
 273.
       
 274.
              for (e_it=mesh_.edges_begin(); e_it!=e_end; ++e_it)
 275.
              {
 276.
                      w  = 0.0;
 277.
       
 278.
                      h0 = mesh_.halfedge_handle(e_it.handle(), 0);
 279.
                      v0 = mesh_.to_vertex_handle(h0);
 280.
                      p0 = mesh_.point(v0);
 281.
       
 282.
                      h1 = mesh_.halfedge_handle(e_it.handle(), 1);
 283.
                      v1 = mesh_.to_vertex_handle(h1);
 284.
                      p1 = mesh_.point(v1);
 285.
       
 286.
                      h2 = mesh_.next_halfedge_handle(h0);
 287.
                      p2 = mesh_.point(mesh_.to_vertex_handle(h2));
 288.
                      d0 = (p0 - p2).normalize();
 289.
                      d1 = (p1 - p2).normalize();
 290.
                      w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, (d0|d1)))));
 291.
       
 292.
                      h2 = mesh_.next_halfedge_handle(h1);
 293.
                      p2 = mesh_.point(mesh_.to_vertex_handle(h2));
 294.
                      d0 = (p0 - p2).normalize();
 295.
                      d1 = (p1 - p2).normalize();
 296.
                      w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, (d0|d1)))));
 297.
       
 298.
                      w = std::max(0.0f, w);
 299.
                      mesh_.property(eweight_,e_it) = w;
 300.
              }
 301.
      }
 302.
       
 303.
      //-----------------------------------------------------------------------------
 304.
       
 305.
      void HarmonicMapViewer::
 306.
      calc_harmonic_parameterization()
 307.
      {
 308.
        // ------------- IMPLEMENT HERE ---------
 309.
              // TASK 5.2 harmonic map computation:
 310.
              // Search and order boundary vertices
 311.
        // Compute boundary parameters
 312.
        // Solve the linear system for internal vertex parameters using solve_linear_system()
 313.
        // Store the parameters in the vparam_h_ vertex property
 314.
              // ------------- IMPLEMENT HERE ---------
 315.
       
 316.
        // Init weights
 317.
        calc_weights();
 318.
       
 319.
        // Calculate parameterization based on eweights_ edge properties
 320.
        gmmVector x,y;
 321.
        calc_parameterization(x,y);
 322.
       
 323.
        // Save to vertices
 324.
        for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){
 325.
          unsigned int i = mesh_.property(vparam_index_,v_it.handle());
 326.
          Vec2f h(x[i],y[i]);
 327.
          mesh_.property(vparam_h_,v_it) = h;
 328.
        }
 329.
      }
 330.
       
 331.
      //-----------------------------------------------------------------------------
 332.
       
 333.
      void
 334.
      HarmonicMapViewer::ComputeTextureCoordinates(int iTextureWidth, int iTextureHeight, int iRepeats, ParameterizationMode imode)
 335.
      {
 336.
        // ------------- IMPLEMENT HERE ---------
 337.
              // TASK 5.3 Compute texture coordinates for textured mesh
 338.
        // rendering using a texture image of dimension iTextureWidth*iTextureHeights and iRepeats repeats
 339.
        // and store them in a mesh property.
 340.
        // If imode is equals to Uniform, compute the texture coordinates using the
 341.
        // parameters stored in vparam_u_ and store the result in texcoord_u_.
 342.
        // If imode is equals to Harmonic, compute the texture coordinates using the
 343.
        // parameters stored in vparam_h_ and store the result in texcoord_h_.
 344.
              // ------------- IMPLEMENT HERE ---------
 345.
       
 346.
        for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){
 347.
       
 348.
          Vec2f uv;
 349.
          if(imode == Uniform)
 350.
            uv = mesh_.property(vparam_u_,v_it);
 351.
          else if(imode == Harmonic)
 352.
            uv = mesh_.property(vparam_h_,v_it);
 353.
          else
 354.
            return;
 355.
       
 356.
          Vec2f uv_n = uv / 2 + Vec2f(.5,.5);
 357.
          uv_n *= iRepeats;
 358.
       
 359.
          if(imode == Uniform)
 360.
            mesh_.property(texcoord_u_,v_it) = uv_n;
 361.
          else if(imode == Harmonic)
 362.
            mesh_.property(texcoord_h_,v_it) = uv_n;
 363.
          else
 364.
            return;
 365.
             
 366.
        }
 367.
      }
 368.
       
 369.
      double HarmonicMapViewer::
 370.
      calc_area(Mesh::FaceHandle fh){
 371.
       
 372.
        // Calculate face area using Heron's formula
 373.
        std::vector<double> lengths;
 374.
        for(Mesh::FaceHalfedgeIter fh_it = mesh_.fh_iter(fh); fh_it; ++fh_it){
 375.
          lengths.push_back(mesh_.calc_edge_length(fh_it.handle()));
 376.
        }
 377.
       
 378.
        assert(lengths.size() == 3);
 379.
       
 380.
        double sp = (lengths[0] + lengths[1] + lengths[2])/2;
 381.
        return sqrt(sp*(sp-lengths[0])*(sp-lengths[1])*(sp-lengths[2]));
 382.
      }
 383.
       
 384.
      double HarmonicMapViewer::
 385.
      calc_area(Vec2f p_a, Vec2f p_b, Vec2f p_c){
 386.
       
 387.
        Vec2f da = p_b - p_a;
 388.
        Vec2f db = p_a - p_c;
 389.
        Vec3f a(da[0],da[1],0);
 390.
        Vec3f b(db[0],db[1],0);
 391.
       
 392.
        return (a % b).length() / 2;
 393.
      }
 394.
       
 395.
      //-----------------------------------------------------------------------------
 396.
       
 397.
      void HarmonicMapViewer::
 398.
      calc_distortion(ParameterizationMode imode)
 399.
      {
 400.
        float angle_distortion=0., area_distortion=0.;
 401.
       
 402.
        // ------------- IMPLEMENT HERE ---------
 403.
              // TASK 5.4 Compute distortion of triangle areas and angles
 404.
        // and print it in the output window.
 405.
        // If imode is equals to Uniform, uniform map distortion has to be
 406.
        // computed.
 407.
        // If imode is equals to Harmonic, harmonic map distortion has to be
 408.
        // computed.
 409.
              // ------------- IMPLEMENT HERE ---------
 410.
       
 411.
       
 412.
        // Calculate total mesh areas
 413.
        Mesh::Scalar total_area_3d = 0;
 414.
        Mesh::Scalar total_area_2d = 0;
 415.
        Mesh::Scalar d1 = 0;
 416.
        for(Mesh::FaceIter f_it=mesh_.faces_begin();f_it!=mesh_.faces_end();++f_it){
 417.
       
 418.
          // Get vertices
 419.
          Mesh::FaceVertexIter fv_it = mesh_.fv_iter(f_it.handle());
 420.
          Mesh::VertexHandle va,vb,vc;
 421.
          va = fv_it.handle();
 422.
          ++fv_it;
 423.
          vb = fv_it.handle();
 424.
          ++fv_it;
 425.
          vc = fv_it.handle();
 426.
       
 427.
          total_area_3d += calc_area(f_it.handle());
 428.
          if(imode == Uniform){
 429.
            total_area_2d += calc_area(mesh_.property(vparam_u_,va),
 430.
                                       mesh_.property(vparam_u_,vb),
 431.
                                       mesh_.property(vparam_u_,vc));
 432.
          } else if(imode == Harmonic){
 433.
            total_area_2d += calc_area(mesh_.property(vparam_h_,va),
 434.
                                       mesh_.property(vparam_h_,vb),
 435.
                                       mesh_.property(vparam_h_,vc));
 436.
          }
 437.
        }
 438.
       
 439.
        // Calculate distortions
 440.
        for(Mesh::FaceIter f_it=mesh_.faces_begin();f_it!=mesh_.faces_end();++f_it){
 441.
       
 442.
          // Get  vertices
 443.
          Mesh::FaceVertexIter fv_it = mesh_.fv_iter(f_it.handle());
 444.
          Mesh::VertexHandle va,vb,vc;
 445.
          va = fv_it.handle();
 446.
          ++fv_it;
 447.
          vb = fv_it.handle();
 448.
          ++fv_it;
 449.
          vc = fv_it.handle();
 450.
       
 451.
          // Calculate area distortion
 452.
          double diff = calc_area(f_it.handle()) / total_area_3d;
 453.
          if(imode == Uniform){
 454.
            diff -= calc_area(mesh_.property(vparam_u_,va),
 455.
                              mesh_.property(vparam_u_,vb),
 456.
                              mesh_.property(vparam_u_,vc)) / total_area_2d;
 457.
          } else if(imode == Harmonic){
 458.
            diff -= calc_area(mesh_.property(vparam_h_,va),
 459.
                              mesh_.property(vparam_h_,vb),
 460.
                              mesh_.property(vparam_h_,vc)) / total_area_2d;
 461.
          }
 462.
          area_distortion += diff*diff;
 463.
       
 464.
          // Calculate angle distortion
 465.
          Vec3f a = (mesh_.point(vc) - mesh_.point(vb)).normalize();
 466.
          Vec3f b = (mesh_.point(va) - mesh_.point(vc)).normalize();
 467.
          Vec3f c = (mesh_.point(va) - mesh_.point(vb)).normalize();
 468.
       
 469.
          double PI = atan(1.)*4;
 470.
          double gamma_3d = acos(std::min(.99f,std::max(-.99f, (-a)|b)));
 471.
          double beta_3d = acos(std::min(.99f,std::max(-.99f, a|c)));
 472.
          double alpha_3d = PI - gamma_3d - beta_3d;
 473.
       
 474.
          Vec2f a2,b2,c2;
 475.
          if(imode == Uniform){
 476.
            a2 = (mesh_.property(vparam_u_, vc) - mesh_.property(vparam_u_, vb)).normalize();
 477.
            b2 = (mesh_.property(vparam_u_, va) - mesh_.property(vparam_u_, vc)).normalize();
 478.
            c2 = (mesh_.property(vparam_u_, va) - mesh_.property(vparam_u_, vb)).normalize();
 479.
          } else if(imode == Harmonic){
 480.
            a2 = (mesh_.property(vparam_h_, vc) - mesh_.property(vparam_h_, vb)).normalize();
 481.
            b2 = (mesh_.property(vparam_h_, va) - mesh_.property(vparam_h_, vc)).normalize();
 482.
            c2 = (mesh_.property(vparam_h_, va) - mesh_.property(vparam_h_, vb)).normalize();
 483.
          }
 484.
         
 485.
       
 486.
          double gamma_2d = acos(std::min(.99f,std::max(-.99f, (-a2)|b2)));
 487.
          double beta_2d = acos(std::min(.99f,std::max(-.99f, a2|c2)));
 488.
          double alpha_2d = PI - gamma_2d - beta_2d;
 489.
       
 490.
          double d_gamma = gamma_3d - gamma_2d;
 491.
          double d_beta = beta_3d - beta_2d;
 492.
          double d_alpha = alpha_3d - alpha_2d;
 493.
          angle_distortion += d_gamma*d_gamma;
 494.
          angle_distortion += d_beta*d_beta;
 495.
          angle_distortion += d_alpha*d_alpha;
 496.
       
 497.
        }
 498.
       
 499.
        cout << "Parameterization distortion: " <<endl;
 500.
        cout << (imode==Uniform? " * Uniform map: ": " * Harmonic map: ") <<endl;
 501.
        cout << " ---- Angle distortion: " << angle_distortion << " Area distortion: " << area_distortion << endl;
 502.
      }
 503.
       
 504.
      //-----------------------------------------------------------------------------
 505.
       
 506.
       
 507.
      static bool _ParameterizationComputed_u = false, _ParameterizationComputed_h = false;
 508.
      static bool _BoundingBox2DComputed = false;
 509.
       
 510.
      void
 511.
      HarmonicMapViewer::
 512.
      draw(const std::string& _draw_mode)
 513.
      {
 514.
        if (indices_.empty())
 515.
        {
 516.
          MeshViewer::draw(_draw_mode);
 517.
          return;
 518.
        }
 519.
        if (_draw_mode == "UV Domain")
 520.
        {
 521.
          if (_ParameterizationMode_!=NoParameterization)
 522.
          {
 523.
            OpenMesh::VPropHandleT<OpenMesh::Vec2f> & TexCoordHandle = (_ParameterizationMode_ == Uniform? vparam_u_: vparam_h_);
 524.
            Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
 525.
            for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
 526.
            {
 527.
              OpenMesh::Vec2f UVCoord = mesh_.property(TexCoordHandle,v_it);
 528.
              mesh_.set_point(v_it, Mesh::Point(-UVCoord[0], UVCoord[1], 0.));
 529.
            }
 530.
            mesh_.update_normals();
 531.
       
 532.
            if (!_BoundingBox2DComputed)
 533.
            {
 534.
              Mesh::ConstVertexIter  v_it(mesh_.vertices_begin()),
 535.
                v_end(mesh_.vertices_end());
 536.
              _bbMin2D = _bbMax2D = mesh_.point(v_it);
 537.
              for (; v_it!=v_end; ++v_it)
 538.
              {
 539.
                _bbMin2D.minimize(mesh_.point(v_it));
 540.
                _bbMax2D.maximize(mesh_.point(v_it));
 541.
              }
 542.
            }
 543.
       
 544.
            set_scene( (Vec3f)(_bbMin2D + _bbMax2D)*0.5, 0.5*(_bbMin2D - _bbMax2D).norm());
 545.
       
 546.
            MeshViewer::draw("Wireframe");
 547.
          }
 548.
        }
 549.
        else
 550.
        {
 551.
          Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
 552.
          for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
 553.
            mesh_.set_point(v_it, mesh_.property(vpos_,v_it));
 554.
          mesh_.update_normals();
 555.
       
 556.
          set_scene( (Vec3f)(_bbMin3D + _bbMax3D)*0.5, 0.5*(_bbMin3D - _bbMax3D).norm());
 557.
       
 558.
          if (_draw_mode == "Textured mesh" && _ParameterizationMode_!=NoParameterization)
 559.
          {
 560.
            float *& TextureCoord = (_ParameterizationMode_ == Uniform? _TextureCoordinates_u: _TextureCoordinates_h);
 561.
            OpenMesh::VPropHandleT<OpenMesh::Vec2f> & TexCoordHandle = (_ParameterizationMode_ == Uniform? texcoord_u_: texcoord_h_);
 562.
            if (!TextureCoord)
 563.
            {
 564.
              int nvertices = mesh_.n_vertices();
 565.
              TextureCoord = new float[2*nvertices];
 566.
       
 567.
              Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
 568.
              int index=0;
 569.
              for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
 570.
              {
 571.
                OpenMesh::Vec2f UVParam = mesh_.property(TexCoordHandle,v_it);
 572.
                TextureCoord[index++] = UVParam[0];
 573.
                TextureCoord[index++] = UVParam[1];
 574.
              }
 575.
       
 576.
            }
 577.
            glEnable( GL_TEXTURE_2D );
 578.
       
 579.
            glEnable(GL_LIGHTING);
 580.
            glShadeModel(GL_SMOOTH);
 581.
       
 582.
            glEnableClientState(GL_VERTEX_ARRAY);
 583.
            glEnableClientState(GL_NORMAL_ARRAY);
 584.
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
 585.
            GL::glVertexPointer(mesh_.points());
 586.
            GL::glNormalPointer(mesh_.vertex_normals());
 587.
            GL::glTexCoordPointer(2, GL_FLOAT, 0, TextureCoord);
 588.
       
 589.
            glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);
 590.
       
 591.
            glDisableClientState(GL_VERTEX_ARRAY);
 592.
            glDisableClientState(GL_NORMAL_ARRAY);
 593.
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
 594.
       
 595.
            glDisable( GL_TEXTURE_2D );
 596.
          }
 597.
          else MeshViewer::draw(_draw_mode);
 598.
        }
 599.
      }
 600.
       
 601.
       
 602.
      //-----------------------------------------------------------------------------
 603.
       
 604.
      void
 605.
      HarmonicMapViewer::
 606.
      keyboard(int key, int x, int y)
 607.
      {
 608.
        switch (toupper(key))
 609.
        {
 610.
        case 'U':
 611.
          {
 612.
            _ParameterizationMode_ = Uniform;
 613.
            if (!_ParameterizationComputed_u)
 614.
            {
 615.
              calc_uniform_parameterization();
 616.
              ComputeTextureCoordinates(256, 256, _Repeat, _ParameterizationMode_);
 617.
              calc_distortion(_ParameterizationMode_);
 618.
       
 619.
              _ParameterizationComputed_u = true;      
 620.
            }
 621.
            glutPostRedisplay();
 622.
            break;
 623.
          }
 624.
        case 'H':
 625.
          {
 626.
            _ParameterizationMode_ = Harmonic;
 627.
            if (!_ParameterizationComputed_h)
 628.
            {
 629.
              calc_harmonic_parameterization();
 630.
              ComputeTextureCoordinates(256, 256, _Repeat, _ParameterizationMode_);
 631.
              calc_distortion(_ParameterizationMode_);
 632.
       
 633.
              _ParameterizationComputed_h = true;      
 634.
            }
 635.
            glutPostRedisplay();
 636.
            break;
 637.
          }
 638.
        default:
 639.
          {
 640.
            MeshViewer::keyboard(key, x, y);
 641.
            break;
 642.
          }
 643.
              }
 644.
      }
 645.
      //-----------------------------------------------------------------------------
 646.
       
 647.
      //=============================================================================

create new paste | create new version of this paste RAW Paste Data
//============================================================================= // // CLASS HarmonicMapViewer - IMPLEMENTATION // //============================================================================= //== INCLUDES ================================================================= #include "HarmonicMapViewer.hh" //== IMPLEMENTATION ========================================================== HarmonicMapViewer:: HarmonicMapViewer(const char* _title, int _width, int _height, int iRepeat) : MeshViewer(_title, _width, _height) { _Repeat = iRepeat; mesh_.request_vertex_colors(); mesh_.add_property(vpos_); mesh_.add_property(vparam_u_); mesh_.add_property(vparam_h_); mesh_.add_property(texcoord_u_); mesh_.add_property(texcoord_h_); mesh_.add_property(vparam_index_); mesh_.add_property(eweight_); add_draw_mode("UV Domain"); add_draw_mode("Textured mesh"); _TextureCoordinates_u=NULL; _TextureCoordinates_h=NULL; _ParameterizationMode_ = NoParameterization; init(); } HarmonicMapViewer:: ~HarmonicMapViewer() { if (glIsTexture(textureID_)) glDeleteTextures( 1, &textureID_); if (_TextureCoordinates_u) delete [] _TextureCoordinates_u; if (_TextureCoordinates_h) delete [] _TextureCoordinates_h; } void HarmonicMapViewer:: init() { // base class first MeshViewer::init(); // generate checkerboard-like image GLubyte tex[256][256][3]; int index=0; for (int x=0; x<256; ++x) { for (int y=0; y<256; ++y) { if ((x<128&&y<128) ||(x>128&&y>128)) { tex[x][y][0] = 0; tex[x][y][1] = 255; tex[x][y][2] = 0; } else { tex[x][y][0] = 255; tex[x][y][1] = 255; tex[x][y][2] = 255; } } } // generate texture glPixelStorei(GL_UNPACK_ALIGNMENT, 1); if (!glIsTexture(textureID_)) glGenTextures(1, &textureID_); glBindTexture(GL_TEXTURE_2D, textureID_); // copy texture to GL glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); glTexImage2D(GL_TEXTURE_2D, 0, 3, 256, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, tex); } bool HarmonicMapViewer:: open_mesh(const char* _meshfilename) { // load mesh if (MeshViewer::open_mesh(_meshfilename)) { // store vertex initial positions and 3D mesh bounding box Mesh::VertexIter v_it=mesh_.vertices_begin(), v_end(mesh_.vertices_end()); _bbMin3D = _bbMax3D = mesh_.point(v_it); for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it) { mesh_.property(vpos_,v_it) = mesh_.point(v_it); _bbMin3D.minimize(mesh_.point(v_it)); _bbMax3D.maximize(mesh_.point(v_it)); } return true; } return false; } //----------------------------------------------------------------------------- void HarmonicMapViewer::solve_linear_system( gmmMatrix& _M, gmmVector& _b, gmmVector& _x) { unsigned int N = _b.size(); _x.resize(N); std::vector< size_t > ipvt(N); gmm::lu_factor( _M, ipvt ); gmm::lu_solve( _M, ipvt, _x, _b ); } void HarmonicMapViewer:: calc_parameterization(gmmVector& x, gmmVector& y){ // Init data size_t n = mesh_.n_vertices(); gmmMatrix W(n,n); gmmVector bx(n),by(n); Mesh::Scalar PI = atan(1.)*4; // Find a border half edge HalfedgeHandle h_init_boundary; bool found_boundary = false; for(Mesh::HalfedgeIter h_it=mesh_.halfedges_begin(); h_it!=mesh_.halfedges_end(); ++h_it){ if(mesh_.is_boundary(h_it)){ h_init_boundary = h_it.handle(); found_boundary = true; break; } } assert(found_boundary); // Calculate total boundary length Mesh::HalfedgeHandle h_cur = h_init_boundary; Mesh::Scalar boundary_length = 0; unsigned int n_boundary = 0; do { boundary_length += mesh_.calc_edge_length(h_cur); h_cur = mesh_.next_halfedge_handle(h_cur); n_boundary++; } while(h_cur != h_init_boundary); // Iterate over boundary h_cur = h_init_boundary; unsigned int i=0; Mesh::Scalar theta=0; do { // Get vertex, save vertex index in matrix Mesh::VertexHandle vh = mesh_.to_vertex_handle(h_cur); mesh_.property(vparam_index_,vh) = i; W(i,i) = 1; bx[i] = cos(theta); by[i] = sin(theta); // Next boundary halfedge h_cur = mesh_.next_halfedge_handle(h_cur); i++; theta += (mesh_.calc_edge_length(h_cur)/boundary_length)*2*PI; } while(h_cur != h_init_boundary); // Set the indices for all non-boundary vertices for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){ if(!mesh_.is_boundary(v_it)){ mesh_.property(vparam_index_,v_it.handle()) = i++; } } // Set matrix values for non-boundary vertices for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){ if(!mesh_.is_boundary(v_it)){ unsigned int i = mesh_.property(vparam_index_,v_it); Mesh::Scalar total_weight = 0; // All outgoing half edges for(Mesh::VertexOHalfedgeIter voh_it=mesh_.voh_iter(v_it.handle()); voh_it;++voh_it){ Mesh::VertexHandle vh_to = mesh_.to_vertex_handle(voh_it.handle()); unsigned int j = mesh_.property(vparam_index_,vh_to); Mesh::EdgeHandle eh = mesh_.edge_handle(voh_it.handle()); Mesh::Scalar weight = mesh_.property(eweight_,eh) / 2; // XXX: Notice division! W(i,j) = weight; total_weight += weight; } W(i,i) = -total_weight; } } //for(unsigned int i=0;i<n;i++){ // for(unsigned int j=0;j<n;j++) // cout << W(i,j) << " "; // // cout << endl; //} // Solve linear systems gmmMatrix W_t(W); solve_linear_system(W,by,y); solve_linear_system(W_t,bx,x); } //----------------------------------------------------------------------------- void HarmonicMapViewer:: calc_uniform_parameterization() { // ------------- IMPLEMENT HERE --------- // TASK 5.1 Uniform map computation: // Search and order boundary vertices // Compute boundary parameters // Solve the linear system for internal vertex parameters using solve_linear_system() // Store the parameters in the vparam_u_ vertex property // ------------- IMPLEMENT HERE --------- // Init weights for(Mesh::EdgeIter e_it=mesh_.edges_begin();e_it!=mesh_.edges_end();++e_it){ // Weight will be divided evenly between vertices, hence resulting vertice weights will be 1! mesh_.property(eweight_, e_it) = 2; } // Calculate parameterization based on eweights_ edge properties gmmVector x,y; calc_parameterization(x,y); // Save to vertices for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){ unsigned int i = mesh_.property(vparam_index_,v_it.handle()); Vec2f u(x[i],y[i]); mesh_.property(vparam_u_,v_it) = u; } } void HarmonicMapViewer:: calc_weights() { Mesh::VertexIter v_it, v_end(mesh_.vertices_end()); Mesh::EdgeIter e_it, e_end(mesh_.edges_end()); Mesh::VertexFaceIter vf_it; Mesh::FaceVertexIter fv_it; Mesh::HalfedgeHandle h0, h1, h2; Mesh::VertexHandle v0, v1; Mesh::Point p0, p1, p2, d0, d1; Mesh::Scalar w; for (e_it=mesh_.edges_begin(); e_it!=e_end; ++e_it) { w = 0.0; h0 = mesh_.halfedge_handle(e_it.handle(), 0); v0 = mesh_.to_vertex_handle(h0); p0 = mesh_.point(v0); h1 = mesh_.halfedge_handle(e_it.handle(), 1); v1 = mesh_.to_vertex_handle(h1); p1 = mesh_.point(v1); h2 = mesh_.next_halfedge_handle(h0); p2 = mesh_.point(mesh_.to_vertex_handle(h2)); d0 = (p0 - p2).normalize(); d1 = (p1 - p2).normalize(); w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, (d0|d1))))); h2 = mesh_.next_halfedge_handle(h1); p2 = mesh_.point(mesh_.to_vertex_handle(h2)); d0 = (p0 - p2).normalize(); d1 = (p1 - p2).normalize(); w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, (d0|d1))))); w = std::max(0.0f, w); mesh_.property(eweight_,e_it) = w; } } //----------------------------------------------------------------------------- void HarmonicMapViewer:: calc_harmonic_parameterization() { // ------------- IMPLEMENT HERE --------- // TASK 5.2 harmonic map computation: // Search and order boundary vertices // Compute boundary parameters // Solve the linear system for internal vertex parameters using solve_linear_system() // Store the parameters in the vparam_h_ vertex property // ------------- IMPLEMENT HERE --------- // Init weights calc_weights(); // Calculate parameterization based on eweights_ edge properties gmmVector x,y; calc_parameterization(x,y); // Save to vertices for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){ unsigned int i = mesh_.property(vparam_index_,v_it.handle()); Vec2f h(x[i],y[i]); mesh_.property(vparam_h_,v_it) = h; } } //----------------------------------------------------------------------------- void HarmonicMapViewer::ComputeTextureCoordinates(int iTextureWidth, int iTextureHeight, int iRepeats, ParameterizationMode imode) { // ------------- IMPLEMENT HERE --------- // TASK 5.3 Compute texture coordinates for textured mesh // rendering using a texture image of dimension iTextureWidth*iTextureHeights and iRepeats repeats // and store them in a mesh property. // If imode is equals to Uniform, compute the texture coordinates using the // parameters stored in vparam_u_ and store the result in texcoord_u_. // If imode is equals to Harmonic, compute the texture coordinates using the // parameters stored in vparam_h_ and store the result in texcoord_h_. // ------------- IMPLEMENT HERE --------- for(Mesh::VertexIter v_it=mesh_.vertices_begin();v_it!=mesh_.vertices_end();++v_it){ Vec2f uv; if(imode == Uniform) uv = mesh_.property(vparam_u_,v_it); else if(imode == Harmonic) uv = mesh_.property(vparam_h_,v_it); else return; Vec2f uv_n = uv / 2 + Vec2f(.5,.5); uv_n *= iRepeats; if(imode == Uniform) mesh_.property(texcoord_u_,v_it) = uv_n; else if(imode == Harmonic) mesh_.property(texcoord_h_,v_it) = uv_n; else return; } } double HarmonicMapViewer:: calc_area(Mesh::FaceHandle fh){ // Calculate face area using Heron's formula std::vector<double> lengths; for(Mesh::FaceHalfedgeIter fh_it = mesh_.fh_iter(fh); fh_it; ++fh_it){ lengths.push_back(mesh_.calc_edge_length(fh_it.handle())); } assert(lengths.size() == 3); double sp = (lengths[0] + lengths[1] + lengths[2])/2; return sqrt(sp*(sp-lengths[0])*(sp-lengths[1])*(sp-lengths[2])); } double HarmonicMapViewer:: calc_area(Vec2f p_a, Vec2f p_b, Vec2f p_c){ Vec2f da = p_b - p_a; Vec2f db = p_a - p_c; Vec3f a(da[0],da[1],0); Vec3f b(db[0],db[1],0); return (a % b).length() / 2; } //----------------------------------------------------------------------------- void HarmonicMapViewer:: calc_distortion(ParameterizationMode imode) { float angle_distortion=0., area_distortion=0.; // ------------- IMPLEMENT HERE --------- // TASK 5.4 Compute distortion of triangle areas and angles // and print it in the output window. // If imode is equals to Uniform, uniform map distortion has to be // computed. // If imode is equals to Harmonic, harmonic map distortion has to be // computed. // ------------- IMPLEMENT HERE --------- // Calculate total mesh areas Mesh::Scalar total_area_3d = 0; Mesh::Scalar total_area_2d = 0; Mesh::Scalar d1 = 0; for(Mesh::FaceIter f_it=mesh_.faces_begin();f_it!=mesh_.faces_end();++f_it){ // Get vertices Mesh::FaceVertexIter fv_it = mesh_.fv_iter(f_it.handle()); Mesh::VertexHandle va,vb,vc; va = fv_it.handle(); ++fv_it; vb = fv_it.handle(); ++fv_it; vc = fv_it.handle(); total_area_3d += calc_area(f_it.handle()); if(imode == Uniform){ total_area_2d += calc_area(mesh_.property(vparam_u_,va), mesh_.property(vparam_u_,vb), mesh_.property(vparam_u_,vc)); } else if(imode == Harmonic){ total_area_2d += calc_area(mesh_.property(vparam_h_,va), mesh_.property(vparam_h_,vb), mesh_.property(vparam_h_,vc)); } } // Calculate distortions for(Mesh::FaceIter f_it=mesh_.faces_begin();f_it!=mesh_.faces_end();++f_it){ // Get vertices Mesh::FaceVertexIter fv_it = mesh_.fv_iter(f_it.handle()); Mesh::VertexHandle va,vb,vc; va = fv_it.handle(); ++fv_it; vb = fv_it.handle(); ++fv_it; vc = fv_it.handle(); // Calculate area distortion double diff = calc_area(f_it.handle()) / total_area_3d; if(imode == Uniform){ diff -= calc_area(mesh_.property(vparam_u_,va), mesh_.property(vparam_u_,vb), mesh_.property(vparam_u_,vc)) / total_area_2d; } else if(imode == Harmonic){ diff -= calc_area(mesh_.property(vparam_h_,va), mesh_.property(vparam_h_,vb), mesh_.property(vparam_h_,vc)) / total_area_2d; } area_distortion += diff*diff; // Calculate angle distortion Vec3f a = (mesh_.point(vc) - mesh_.point(vb)).normalize(); Vec3f b = (mesh_.point(va) - mesh_.point(vc)).normalize(); Vec3f c = (mesh_.point(va) - mesh_.point(vb)).normalize(); double PI = atan(1.)*4; double gamma_3d = acos(std::min(.99f,std::max(-.99f, (-a)|b))); double beta_3d = acos(std::min(.99f,std::max(-.99f, a|c))); double alpha_3d = PI - gamma_3d - beta_3d; Vec2f a2,b2,c2; if(imode == Uniform){ a2 = (mesh_.property(vparam_u_, vc) - mesh_.property(vparam_u_, vb)).normalize(); b2 = (mesh_.property(vparam_u_, va) - mesh_.property(vparam_u_, vc)).normalize(); c2 = (mesh_.property(vparam_u_, va) - mesh_.property(vparam_u_, vb)).normalize(); } else if(imode == Harmonic){ a2 = (mesh_.property(vparam_h_, vc) - mesh_.property(vparam_h_, vb)).normalize(); b2 = (mesh_.property(vparam_h_, va) - mesh_.property(vparam_h_, vc)).normalize(); c2 = (mesh_.property(vparam_h_, va) - mesh_.property(vparam_h_, vb)).normalize(); } double gamma_2d = acos(std::min(.99f,std::max(-.99f, (-a2)|b2))); double beta_2d = acos(std::min(.99f,std::max(-.99f, a2|c2))); double alpha_2d = PI - gamma_2d - beta_2d; double d_gamma = gamma_3d - gamma_2d; double d_beta = beta_3d - beta_2d; double d_alpha = alpha_3d - alpha_2d; angle_distortion += d_gamma*d_gamma; angle_distortion += d_beta*d_beta; angle_distortion += d_alpha*d_alpha; } cout << "Parameterization distortion: " <<endl; cout << (imode==Uniform? " * Uniform map: ": " * Harmonic map: ") <<endl; cout << " ---- Angle distortion: " << angle_distortion << " Area distortion: " << area_distortion << endl; } //----------------------------------------------------------------------------- static bool _ParameterizationComputed_u = false, _ParameterizationComputed_h = false; static bool _BoundingBox2DComputed = false; void HarmonicMapViewer:: draw(const std::string& _draw_mode) { if (indices_.empty()) { MeshViewer::draw(_draw_mode); return; } if (_draw_mode == "UV Domain") { if (_ParameterizationMode_!=NoParameterization) { OpenMesh::VPropHandleT<OpenMesh::Vec2f> & TexCoordHandle = (_ParameterizationMode_ == Uniform? vparam_u_: vparam_h_); Mesh::VertexIter v_it, v_end(mesh_.vertices_end()); for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it) { OpenMesh::Vec2f UVCoord = mesh_.property(TexCoordHandle,v_it); mesh_.set_point(v_it, Mesh::Point(-UVCoord[0], UVCoord[1], 0.)); } mesh_.update_normals(); if (!_BoundingBox2DComputed) { Mesh::ConstVertexIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); _bbMin2D = _bbMax2D = mesh_.point(v_it); for (; v_it!=v_end; ++v_it) { _bbMin2D.minimize(mesh_.point(v_it)); _bbMax2D.maximize(mesh_.point(v_it)); } } set_scene( (Vec3f)(_bbMin2D + _bbMax2D)*0.5, 0.5*(_bbMin2D - _bbMax2D).norm()); MeshViewer::draw("Wireframe"); } } else { Mesh::VertexIter v_it, v_end(mesh_.vertices_end()); for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it) mesh_.set_point(v_it, mesh_.property(vpos_,v_it)); mesh_.update_normals(); set_scene( (Vec3f)(_bbMin3D + _bbMax3D)*0.5, 0.5*(_bbMin3D - _bbMax3D).norm()); if (_draw_mode == "Textured mesh" && _ParameterizationMode_!=NoParameterization) { float *& TextureCoord = (_ParameterizationMode_ == Uniform? _TextureCoordinates_u: _TextureCoordinates_h); OpenMesh::VPropHandleT<OpenMesh::Vec2f> & TexCoordHandle = (_ParameterizationMode_ == Uniform? texcoord_u_: texcoord_h_); if (!TextureCoord) { int nvertices = mesh_.n_vertices(); TextureCoord = new float[2*nvertices]; Mesh::VertexIter v_it, v_end(mesh_.vertices_end()); int index=0; for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it) { OpenMesh::Vec2f UVParam = mesh_.property(TexCoordHandle,v_it); TextureCoord[index++] = UVParam[0]; TextureCoord[index++] = UVParam[1]; } } glEnable( GL_TEXTURE_2D ); glEnable(GL_LIGHTING); glShadeModel(GL_SMOOTH); glEnableClientState(GL_VERTEX_ARRAY); glEnableClientState(GL_NORMAL_ARRAY); glEnableClientState(GL_TEXTURE_COORD_ARRAY); GL::glVertexPointer(mesh_.points()); GL::glNormalPointer(mesh_.vertex_normals()); GL::glTexCoordPointer(2, GL_FLOAT, 0, TextureCoord); glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]); glDisableClientState(GL_VERTEX_ARRAY); glDisableClientState(GL_NORMAL_ARRAY); glDisableClientState(GL_TEXTURE_COORD_ARRAY); glDisable( GL_TEXTURE_2D ); } else MeshViewer::draw(_draw_mode); } } //----------------------------------------------------------------------------- void HarmonicMapViewer:: keyboard(int key, int x, int y) { switch (toupper(key)) { case 'U': { _ParameterizationMode_ = Uniform; if (!_ParameterizationComputed_u) { calc_uniform_parameterization(); ComputeTextureCoordinates(256, 256, _Repeat, _ParameterizationMode_); calc_distortion(_ParameterizationMode_); _ParameterizationComputed_u = true; } glutPostRedisplay(); break; } case 'H': { _ParameterizationMode_ = Harmonic; if (!_ParameterizationComputed_h) { calc_harmonic_parameterization(); ComputeTextureCoordinates(256, 256, _Repeat, _ParameterizationMode_); calc_distortion(_ParameterizationMode_); _ParameterizationComputed_h = true; } glutPostRedisplay(); break; } default: { MeshViewer::keyboard(key, x, y); break; } } } //----------------------------------------------------------------------------- //=============================================================================
Pastebin.com Tools & Applications
Windows Desktop
Firefox
Chrome
iPhone & iPad
HP WebOS
Android
Mac OS
Opera
pastebin.com
create new paste | API | trends | faq | tools | privacy | contact | stats | in the news | go pro
pastebin domains center | whats new in pastebin v3? | public todo list | advertise on pastebin
Our Sites: GoogleCounter | W3Wiki | Hostlogr | Tinysubs | Server Time: 0.01224
