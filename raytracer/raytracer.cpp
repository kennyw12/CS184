#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <time.h>
#include <math.h>
#include <string.h>

struct Point
{
    int x, y, z;
    Point();
    Point(int x1, int y1, int z1) {
        x = x1;
        y = y1;
        z = z1;
    }
};

struct Vector
{
    int x, y, z;
    float mag;
    Vector(int x1, int y1, int z1) {
        x = x1;
        y = y1;
        z = z1;
        mag = pow(x1*x1 + y1*y1 + z1*z1, 1/2);
    }
    Vector(Point* a, Point* b) {
        x = a->x - b->x;
        y = a->y - b->y;
        z = a->z - b->z;
        mag = pow(x*x + y*y + z*z, 1/2);
    }
};


Vector* normalize(Vector* v) {
  return new Vector(v->x/v->mag, v->y/v->mag, v->z/v->mag);
}

struct Ray
{
    Vector* direction;
    Point* origin;
    Ray(Point* p, Vector* v) {
        origin = p;
        direction = normalize(v);
    }
};

int dotProduct(Vector* a, Vector* b) {
  return a->x * b->x + a->y * b->y + a->z * b->z;
}

Vector* subtract(Vector* a, Vector* b) {
  return new Vector(a->x - b->x, a->y - b->y, a->z - b->z);
}

Vector* add(Vector* a, Vector* b) {
  return new Vector(a->x + b->x, a->y + b->y, a->z + b->z);
}

Vector* scale(Vector* a, int b) {
  return new Vector(a->x * b, a->y * b, a->z * b);
}

Vector* reflect(Vector* l, Vector* n) {
  float dot = dotProduct(l, n);
  return add(scale(scale(n, dot), 2), scale(l, -1));
}


class Shape {
    protected:
        Point* center;
    public:
        void set_center(int a, int b, int c) {
            center = new Point(a, b, c);
        }
        virtual float get_color(Ray* ray) = 0;
        virtual int get_radius() = 0;
};

class Sphere: public Shape {
    int radius;
    public:
        Sphere(int x1, int y1, int z1, int r) {
            //printf("Set Radius %d\n", r);
            set_center(x1, y1, z1);
            radius = r;
        }
        float get_color(Ray* ray) {
            //printf("Getting Color %d\n", radius);
            Vector* diff = new Vector(ray->origin, center);
            float det = pow(dotProduct(ray->direction, diff), 2) - dotProduct(diff,diff) + (radius*radius);
            printf("det %f\n", det);
            if (det < 0) {
                return 0; //miss
            } else if (det == 0) {
                //color computation for a Vector
                return 1;
            } else if (det > 0) {
                //2 intersections
                return 1;
            }
        }
        int get_radius() {
            //printf("Getting Radius\n");
            return radius;
        }
};


struct Scene
{
    Point* ll;
    Point* lr;
    Point* ul;
    Point* ur;
    Point* eye;
    int width, height, depth;
    float* rsamples;
    float* gsamples;
    float* bsamples;
    std::vector<Shape*> shapes;
    Scene(Point* ll1, Point* lr1, Point* ul1, Point* ur1, Point* eye1, int w, int h) {
        ll = ll1;
        lr = lr1;
        ul = ul1;
        ur = ur1;
        eye = eye1;
        width = w;
        height = h;
        depth = (ll->z + lr->z + ul->z + ur->z)/4;
        rsamples = new float[w*h];
        gsamples = new float[w*h];
        bsamples = new float[w*h];
    }
    void shoot_ray(Ray* ray, int i, int j) {
        Shape* s = shapes.at(0);
        float color = s->get_color(ray);
        rsamples[i + j*height] = color;
    }
    void add_shape(Shape* s) {
        shapes.push_back(s);
    }
};

struct Sampler
{
    Scene* my_scene;
    Sampler(Scene* scene) {
        my_scene = scene;
    }
    void sample() {
        int width = my_scene->width;
        int height = my_scene->height;
        int depth = my_scene->depth;
        for(int i=0;i<width;i++) {
            for(int j=0;j<height;j++) {
                //printf("i %d j %d\n", i, j);
                Point* pos = new Point(i - width/2, j - height/2, depth);
                //printf("About to make Vector\n");
                Ray* ray = new Ray(my_scene->eye, new Vector(pos, my_scene->eye));
                //printf("About to Generate\n");
                my_scene->shoot_ray(ray, i, j);
            }
        }
    }
};

struct Film
{
    Sampler* my_sampler;
    Film(Sampler* sampler) {
        my_sampler = sampler;
    }
    void writeImage() {
        int width = my_sampler->my_scene->width;
        int height = my_sampler->my_scene->height;
        my_sampler->sample();
        FILE *file = fopen("image.ppm", "wb");
        fprintf(file, "P3\n%d %d\n%d\n", width, height, 255);

        int count = 0;
        int x, y;
        // Write image data
        for (y=0 ; y<height ; y++) {
            for (x=0 ; x<width ; x++) {
               count += 1;
               int r = my_sampler->my_scene->rsamples[y*width + x] * 255;
               int g = my_sampler->my_scene->gsamples[y*width + x] * 255;
               int b = my_sampler->my_scene->bsamples[y*width + x] * 255;
               fprintf(file, " %d %d %d    ", r, g, b);
               if (count >= 5) {
                   fprintf(file, "\n");
                   count = 0;
               }
            }
        }
        fclose(file);
    }
};


int main(int argc, char *argv[]) {
    int i = 1;
    int w = 100;
    int h = 100;
    Point* ll = new Point(-1,-1,-1);
    Point* lr = new Point(1,-1,-1);
    Point* ul = new Point(-1,1,-1);
    Point* ur = new Point(1,1,-1);
    Point* eye = new Point(0,0,0);
    printf("Start\n");
    //printf("%d %d %d\n", argc, i, i<argc);
    while (i < argc) {
        if (!(strcmp(argv[i],"-eye"))) {
            eye = new Point(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
            i += 4;
        }
        else if (!(strcmp(argv[i],"-ll"))) {
            ll = new Point(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
            i += 4;
        }
        else if (!(strcmp(argv[i],"-lr"))) {
            lr = new Point(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
            i += 4;
        }
        else if (!(strcmp(argv[i],"-ul"))) {
            ul = new Point(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
            i += 4;
        }
        else if (!(strcmp(argv[i],"-ur"))) {
            ur = new Point(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
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
    Scene* scene = new Scene(ll, lr, ul, ur, eye, w, h);
    scene->add_shape((Shape*) new Sphere(0,0,-2,1));
    Sampler* sampler = new Sampler(scene);
    Film* my_film = new Film(sampler);
    my_film->writeImage();
    return 0;
}