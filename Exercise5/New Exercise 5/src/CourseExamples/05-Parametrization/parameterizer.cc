                                                                          //=============================================================================
#include "HarmonicMapViewer.hh"



int main(int argc, char **argv)
{
  if (argc < 3) 
	{
		std::cerr << "Usage: \n" 
			<< argv[0] << " InputMesh.off Repeats \n\n";
		exit(1);
	}

  glutInit(&argc, argv);

  HarmonicMapViewer window("Harmonic Map", 512, 512, atoi(argv[2]));

  window.open_mesh(argv[1]);

  glutMainLoop();
}
