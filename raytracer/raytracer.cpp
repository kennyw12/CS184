#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <time.h>
#include <math.h>
#include <string.h>
#include "raytracer.h"

void vector::vector(int ix, int iy, int iz) {
	x = ix;
	y = iy;
	z = iz;
	magnitude = pow((x*x + y*y + z*z), 0.5);
}	

rgb::rgb(float ir, float ig, float ib) {
    r = ir;
	g = ig;
	b = ib;
}

scene::scene(vector ieye, vector ll, vector lr, vector ul, vector ur, int w, int h) {
	eye = ieye;
	lower_left = ll;
	lower_right = lr;
	upper_left = ul;
	upper_right = ur;
	width = w;
	height = h;
	depth = (ll.z + lr.z + ul.z + ur.z)/4;
}

sampler::sampler(scene s1) {
	my_scene = s1;
}
rgb sampler::generate_sample(vector ray) {
	return rgb(0,0,0);
}

film::film(scene s1, sampler sample) {
	my_scene = s1;
	my_sampler = sample;
	rgb s[my_scene.height][my_scene.width];
	samples = s;
}
void film::sample() {
	for(int i=0;i<my_scene.width;i++) {
		for(int j=0;j<my_scene.height;j++) {
			vector ray = make_ray(my_scene.eye, (i, j, my_scene.depth));
			samples[j][i] = my_scene.generate_sample(ray);
		}
	}	
}
void film::write_to_file() {
	FILE file = fopen("image.ppm", "wb");
	fprintf(fp, "P6\n%d %d\n255\n", my_scene.width, my_scene.height);
	for (int i = 0; i < my_scene.height; i+=1) {
	    for (int j = 0; j < my_scene.width; j+=1) {
	    	static unsigned char color[3];
	    	rgb pix = samples[i][j];
		    color[0] = pix * 256;  /* red */
		    color[1] = pix * 256;  /* green */
		    color[2] = pix * 256;  /* blue */
		    fwrite(color, 1, 3, file);
	    }
	}
	fclose(file);
}

vector make_ray(vector v1, vector v2) {
	return vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

void main(int argc, char *argv[]) {
	int i = 0;
	vector eye = new vector(0,0,0);
	vector ll = new vector(-1,-1,-1);
	vector lr = new vector(-1,1,-1);
	vector ul = new vector(1,-1,-1);
	vector ur = new vector(1,1,-1);
	int w = 100;
	int h = 100;
	while (i < argc) {
		if (!(strcmp(argv[i],"-eye"))) {
			eye = new vector(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
			i += 4;
		}
		else if (!(strcmp(argv[i],"-ll"))) {
			ll = new vector(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
			i += 4;
		}
		else if (!(strcmp(argv[i],"-lr"))) {
			lr = new vector(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
			i += 4;
		}
		else if (!(strcmp(argv[i],"-ul"))) {
			ul = new vector(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
			i += 4;
		}
		else if (!(strcmp(argv[i],"-ur"))) {
			ur = new vector(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
			i += 4;
		}
		else if (!(strcmp(argv[i],"-w"))) {
			w = atof(argv[i+1]);
			i += 2;
		}
		else if (!(strcmp(argv[i],"-h"))) {
			h = atof(argv[i+1]);
			i += 2;
		}
	}
	scene my_scene = new scene(eye, ll, lr, ul, ur, w, h);
	sampler my_sampler = new sampler(my_scene);
	film my_film = new film(my_scene. my_sampler);
	film.sample();
	film.write_to_file();
}