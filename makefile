minimalsurfaces:
	g++ -std=gnu++11 -o minimal_surfaces src\minimal_surfaces.cpp src\mesh.cpp src\example_surfaces.cpp -lglfw3 -lgdi32 -lglu32 -lglew32  -lopengl32 -IC:\eigen-eigen\