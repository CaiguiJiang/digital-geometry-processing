Color Coding:
=-=-=-=-=-=-=
The color is selected according to the valence, according to the following table:

Valence		Color		RGB
-------		-----		---
1		Violet 		(255,0,255)
2		White		(255,255,255)
3		Red		(255,0,0)
4		Orange		(255,125,0)
5		Yellow		(255,255,0)
6		Green		(0,255,0)
7		Cyan		(0,255,255)
8		Blue		(0,0,255)
9		Indigo		(125,0,125)
10+		Maroon		(125,0,0)

Implementation Details:
=-=-=-=-=-=-=-=-=-=-=-=
A custom vertex property is defined in ValenceViewer.hh - "vertexValence".
It is of type int, and it is used for storing the calculated valence of each vertex.
The color coding function is implemented using a mapping from int to Mesh::Color -
The key being the vertex's valence, and the value being the color, according to the above table.
The mapping is implemented using an array.
The array is initialized in the ValenceViewer::InitializeColorsArray() method, which is called in the c'tor.