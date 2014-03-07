class vector {
		int x, y, z;
		float magnitude;
	public:
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
	scene(vector ieye, vector ll, vector lr, vector ul, vector ur, int w, int h);
};

class sampler {
	scene my_scene;
public:
	sampler(scene s1);
	rgb generate_sample(vector ray);
};

class film {
	scene my_scene;
	sampler my_sampler;
	rgb* samples;
	const int max_color_value = 255;
public:
	film(scene s1, sampler sample);
	void sample();
	void write_to_file();
};