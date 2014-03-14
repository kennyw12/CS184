class vector {
		int x, y, z;
		float magnitude;
	public
		vector(int ix, int iy, int iz);
};

class rgb {
	float r, g, b;
public:
	rgb(float ir, float ig, float ib);
};

class scene {
	int width, height, depth;
	vector eye, lower_left, lower_right, upper_left, upper_right;
public:
	scene(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int);
};

class sampler {
	scene my_scene;
public:
	sampler(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int);
	rgb generate_sample(vector ray);
};

class film {
	scene my_scene;
	sampler my_sampler;
	rgb* samples;
	const int max_color_value = 255;
public:
	film(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int);
	void sample();
	void write_to_file();
};