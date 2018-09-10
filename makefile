minimalsurfaces:
	g++ -std=gnu++11 -o minimalsurfaces MinimalSurfaces.cpp mesh.cpp example_surfaces.cpp -lfreeglut -lopengl32 -lglu32 -lglfw3 -lglew32 -IC:\eigen-eigen\