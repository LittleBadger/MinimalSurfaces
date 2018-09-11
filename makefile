minimalsurfaces:
	g++ -std=gnu++11 -o minimal_surfaces src\minimal_surfaces.cpp src\mesh.cpp src\example_surfaces.cpp -lfreeglut -lopengl32 -lglu32 -lglfw3 -lglew32 -IC:\eigen-eigen\