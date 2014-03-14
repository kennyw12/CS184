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

using namespace std;

int width, height, depth;
float zero = 0.0;


struct Point
{
    float x, y, z;
    Point();
    Point(float x1, float y1, float z1) {
        x = x1;
        y = y1;
        z = z1;
    }
};

struct Vector
{
    float x, y, z;
    float mag;
    Vector(float x1, float y1, float z1) {
        x = x1;
        y = y1;
        z = z1;
        mag = pow(x1*x1 + y1*y1 + z1*z1, 0.5);
    }
    Vector(Point* a, Point* b) {
        x = a->x - b->x;
        y = a->y - b->y;
        z = a->z - b->z;
        mag = pow(x*x + y*y + z*z, 0.5);
    }
    Vector(Point* a) {
        x = a->x;
        y = a->y;
        z = a->z;
        mag = pow(x*x + y*y + z*z, 0.5);
    }
};

Point* make_point(Vector* v, Point* p) {
    return new Point(v->x - p->x, v->y - p->y, v->z - p->z);
}

Vector* normalize(Vector* v) {
  return new Vector(v->x/v->mag, v->y/v->mag, v->z/v->mag);
}

Point* normalize(Point* v) {
    float mag = pow(v->x*v->x + v->y*v->y + v->z*v->z, 1/2);
    return new Point(v->x/mag, v->y/mag, v->z/mag);
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

float dotProduct(Vector* a, Vector* b) {
  return a->x * b->x + a->y * b->y + a->z * b->z;
}

float dotProduct(Point* a, Point* b) {
  return a->x * b->x + a->y * b->y + a->z * b->z;
}

float dotProduct(Point* a, Vector* b) {
  return a->x * b->x + a->y * b->y + a->z * b->z;
}

float dotProduct(Vector* a, Point* b) {
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

Vector* scale(Vector* a, float b) {
  return new Vector(a->x * b, a->y * b, a->z * b);
}

Point* subtract(Point* a, Point* b) {
  return new Point(a->x - b->x, a->y - b->y, a->z - b->z);
}

Point* add(Point* a, Point* b) {
  return new Point(a->x + b->x, a->y + b->y, a->z + b->z);
}

Point* scale(Point* a, int b) {
  return new Point(a->x * b, a->y * b, a->z * b);
}

Point* scale(Point* a, float b) {
  return new Point(a->x * b, a->y * b, a->z * b);
}

Vector* reflect(Vector* l, Vector* n) {
  float dot = dotProduct(l, n);
  return add(scale(scale(n, dot), 2), scale(l, -1));
}

struct Color {
  float r, g, b;
  Color(float r1=0.0, float g1=0.0, float b1=0.0) {
    r = r1;
    g = g1;
    b = b1;
  }
};

struct Light {
    Point* pos;
    Color* color;
    Light(Point* p, Color* c) {
        pos = new Point(p->x * width/2, p->y * height/2, p->z * depth/2);
        color = c;
    }
};

std::vector<Light*> plights;
std::vector<Light*> dlights;
Vector* viewer = new Vector(0,0,1);

class Shape {
    protected:
        Point* center;
        Color* ka;
        Color* ks;
        Color* kd;
        int p;
    public:
        void set_center(int a, int b, int c) {
            center = new Point(a, b, c);
        }
        void set_ka(Color* c) {
            ka = c;
        }
        void set_ks(Color* c, int p1) {
            ks = c;
            p = p1;
        }
        void set_kd(Color* c) {
            kd = c;
        }
        virtual Color* get_color(Ray* ray) = 0;
        virtual int get_radius() = 0;
};



class Sphere: public Shape {
    int radius;
    public:
        Sphere(int x1, int y1, int z1, int r, Color* ka, Color* kd, Color* ks, int p) {
            //printf("Set Radius %d\n", r);
            set_center(x1 * width/2, y1 * height/2, z1 * depth/2);
            radius = r * std::min(width, height) /2;
            set_kd(kd);
            set_ks(ks, p);
            set_ka(ka);
        }
        Color* get_color(Ray* ray) {
            //printf("Getting Color %d\n", radius);
            Ray* tray = new Ray(add(ray->origin, center), ray->direction);
            //printf("a\n");
            float a = dotProduct(tray->direction, tray->direction);
            //printf("b\n");
            float b = 2 * dotProduct(tray->direction, tray->origin);
            //printf("c\n");
            float c = dotProduct(tray->origin, tray->origin) - (radius * radius);

            float det = b * b - 4 * a * c;
            //printf("det, %f\n", det);            
            if (det < 0) {
                return new Color(0, 0, 0); //Misses sphere
            } else {
                //color computation for a Vector
                //printf("COLOR\n");
                float t = min((-b-pow(det,0.5))/(2*a), (-b+pow(det,0.5))/(2*a));
                //printf("%f t\n", t);
                Point* intersect = make_point(scale(ray->direction, t), ray->origin);
                //printf("%f, %f, %f %f normal\n", intersect->x, intersect->y, intersect->z, t);
                Vector* normal = normalize(new Vector(intersect, center));
                float r = 0.0;
                float g = 0.0;
                float b = 0.0;
                //printf("light calc, %d, %d\n", plights.size(), dlights.size());
                int i;
                for(i=0;i<plights.size();i++) {
                    Light* light = plights.at(i);
                    Vector* ltos = normalize(new Vector(light->pos, intersect));
                    Vector* reflection = normalize(reflect(ltos, normal));
                    //printf("%f %f %f mags\n", ltos->mag, reflection->mag, inter->mag);
                    //printf("%d, %d, %d\n", ray->origin->x, ray->origin->y, ray->origin->z);
                    //printf("%f, %f, %f ray\n", light->pos->x, light->pos->y, light->pos->z);
                    //printf("%f, %f, %f ltos\n", ltos->x, ltos->y, ltos->z);
                    //printf("%f, %f, %f refl\n", reflection->x, reflection->y, reflection->z);
                    //printf("%f, %f, %f normal\n", normal->x, normal->y, normal->z);
                    float diff = max(dotProduct(ltos, normal), zero);
                    float stuff = dotProduct(viewer, reflection);
                    float spec = pow(max(stuff, zero), p);
                    printf("%f, %f %f pl\n", diff, spec, stuff);
                    r = r + ka->r*light->color->r + kd->r*light->color->r*diff + ks->r*light->color->r*spec;
                    g = g + ka->g*light->color->g + kd->g*light->color->g*diff + ks->g*light->color->g*spec;
                    b = b + ka->b*light->color->b + kd->b*light->color->b*diff + ks->b*light->color->b*spec;
                }
                for(i=0;i<dlights.size();i++) {
                    Light* light = dlights.at(i);
                    Vector* negLight = normalize(new Vector(scale(light->pos, -1)));
                    Vector* reflection = normalize(reflect(negLight, normal));
                    float diff = max(dotProduct(negLight, normal), zero);
                    float stuff = dotProduct(viewer, reflection);
                    float spec = pow(max(stuff, zero), p);
                    printf("%f, %f %f dl\n", diff, spec, stuff);
                    r = r + ka->r*light->color->r + kd->r*light->color->r*diff + ks->r*light->color->r*spec;
                    g = g + ka->g*light->color->g + kd->g*light->color->g*diff + ks->g*light->color->g*spec;
                    b = b + ka->b*light->color->b + kd->b*light->color->b*diff + ks->b*light->color->b*spec;
                }
                //printf("%f, %f, %f\n", r, g, b);
                return new Color(r, g, b);
            }
        }
        int get_radius() {
            //printf("Getting Radius\n");
            return radius;
        }
        Ray* transform_ray(Ray* ray) {

        }
};


struct Scene
{
    Point* ll;
    Point* lr;
    Point* ul;
    Point* ur;
    Point* eye;
    float* rsamples;
    float* gsamples;
    float* bsamples;
    std::vector<Shape*> shapes;
    Scene(Point* ll1, Point* lr1, Point* ul1, Point* ur1, Point* eye1, int w, int h) {
        width = w;
        height = h;
        depth = (pow(w*w + h*h, 0.5))*((ll1->z + lr1->z + ul1->z + ur1->z)/4);
        //printf("%d\n", depth);
        ll = new Point(ll1->x * w/2, ll1->y * h/2, ll1->z * depth/2);
        lr = new Point(lr1->x * w/2, lr1->y * h/2, lr1->z * depth/2);
        ul = new Point(ul1->x * w/2, ul1->y * h/2, ul1->z * depth/2);
        ur = new Point(ur1->x * w/2, ur1->y * h/2, ur1->z * depth/2);
        eye = new Point(eye1->x * w/2, eye1->y * h/2, eye1->z * depth/2);
        rsamples = new float[w*h];
        gsamples = new float[w*h];
        bsamples = new float[w*h];
    }
    void shoot_ray(Ray* ray, int i, int j) {
        Shape* s = shapes.at(0);
        Color* color = s->get_color(ray);
        rsamples[i + j*height] = color->r;
        gsamples[i + j*height] = color->g;
        bsamples[i + j*height] = color->b;
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
        for(int i=0;i<width;i++) {
            for(int j=0;j<height;j++) {
                //printf("i %d j %d\n", i, j);
                Point* pos = new Point(i - width/2, (j - height/2)*-1, depth/2);
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
    int w = 400;
    int h = 400;
    Point* ll = new Point(-1,-1,-1);
    Point* lr = new Point(1,-1,-1);
    Point* ul = new Point(-1,1,-1);
    Point* ur = new Point(1,1,-1);
    Point* eye = new Point(0,0,0);
    //printf("Start\n");
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
    Color* ka = new Color(0.1, 0.1, 0);
    Scene* scene = new Scene(ll, lr, ul, ur, eye, w, h);
    scene->add_shape((Shape*) new Sphere(0,0,-2,1, ka, new Color(1, 1, 0), new Color(.8,.8,.8), 16));
    //printf("lights\n");
    plights.push_back(new Light(new Point(2, 2, 0), new Color(.6, .6, .6)));
    //dlights.push_back(new Light(new Point(10, 10, 10), new Color(1, 1, 1)));
    //printf("sampler\n");
    Sampler* sampler = new Sampler(scene);
    Film* my_film = new Film(sampler);
    my_film->writeImage();
    return 0;
}