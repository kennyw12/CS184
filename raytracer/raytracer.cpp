#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <cfloat>

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

Point* add(Point* a, Vector* b) {
  return new Point(a->x + b->x, a->y + b->y, a->z + b->z);
}

Point* scale(Point* a, int b) {
  return new Point(a->x * b, a->y * b, a->z * b);
}

Point* scalez(Point* a, int b) {
  return new Point(a->x, a->y, a->z * b);
}

Vector* scalez(Vector* a, int b) {
  return new Vector(a->x, a->y, a->z * b);
}

Point* scale(Point* a, float b) {
  return new Point(a->x * b, a->y * b, a->z * b);
}

Vector* reflect(Vector* l, Vector* n) {
  float dot = dotProduct(l, n);
  return add(scale(scale(n, dot), 2), scale(l, -1));
}

Vector* crossProduct(Vector* v1, Vector* v2) {
    return new Vector((v1->y*v2->z) - (v1->z*v2->y), (v1->z*v2->x) - (v1->x*v2->z), (v1->x*v2->y) - (v1->y*v2->x));
}

Point* ray_endpoint(Ray* ray, float t) {
    return make_point(scale(ray->direction, t), ray->origin);
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
        pos = new Point(p->x*width/2, p->y*height/2, p->z*depth/3);
        color = c;
    }
};

std::vector<Light*> plights;
std::vector<Light*> dlights;

class Shape {
    public:
        Point* center;
        Color* ka;
        Color* ks;
        Color* kd;
        Color* kr;
        int p;
        int sid;
        int radius;
        void set_center(float a, float b, float c) {
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
        void set_kr(Color* c) {
            kr = c;
        }
        int get_radius() {
            return radius;
        }
};

class Triangle: public Shape {
    public:
        Point* vertices[3];
        Vector* normal;
        Triangle(Point* p1, Point* p2, Point* p3, Color* ka, Color* kd, Color* ks, int p, Color* kr, int sid1) {
            vertices[0] = new Point(p1->x * width/2 * -1, p1->y * height/2 * -1, p1->z * depth/3);
            vertices[1] = new Point(p2->x * width/2 * -1, p2->y * height/2 * -1, p2->z * depth/3);
            vertices[2] = new Point(p3->x * width/2 * -1, p3->y * height/2 * -1, p3->z * depth/3);
            //printf("%f %f %f\n", vertices[0]->x, vertices[0]->y, vertices[0]->z);
            //printf("%f %f %f\n", vertices[1]->x, vertices[1]->y, vertices[1]->z);
            //printf("%f %f %f\n", vertices[2]->x, vertices[2]->y, vertices[2]->z);
            set_center((p1->x+p2->x+p3->x)*width/6,(p1->y+p2->y+p3->y)*height/6,(p1->z+p2->z+p3->z)*depth/3);
            radius = ((float)min(width, height)) * 0.5;
            Vector* v1 = new Vector(p2, p1);
            Vector* v2 = new Vector(p3, p1);
            normal = normalize(crossProduct(v1, v2));
            set_kd(kd);
            set_ks(ks, p);
            set_ka(ka);
            set_kr(kr);
            sid = sid1;
        }
};

class Sphere: public Shape {
    public:
        float radius;
        Sphere(float x1, float y1, float z1, float r, Color* ka, Color* kd, Color* ks, int p, Color* kr, int sid1) {
            //printf("Set Radius %d\n", r);
            set_center(x1 * width/2 * -1, y1 * height/2 * -1, z1 * depth/3);
            printf("%f %f %f\n", center->x, center->y, center->z);
            radius = r * ((float)min(width, height)) * 0.5;
            //printf("%f\n", radius);
            set_kd(kd);
            set_ks(ks, p);
            set_ka(ka);
            set_kr(kr);
            sid = sid1;
        }
};

float time_to_sphere(Ray* ray, Sphere* s) {
    Ray* tray = new Ray(add(ray->origin, s->center), ray->direction);
    //printf("a\n");
    float a = dotProduct(tray->direction, tray->direction);
    //printf("b\n");
    float b = 2 * dotProduct(tray->direction, tray->origin);
    //printf("c\n");
    float c = dotProduct(tray->origin, tray->origin) - (s->radius * s->radius);

    float det = b * b - 4 * a * c;
    //printf("det, %f\n", det);
    if (det<0) {
        return det; //Misses sphere
    } else {
        float t = min((-b-pow(det,0.5))/(2*a), (-b+pow(det,0.5))/(2*a));
        return t;
    }
}
float time_to_triangle(Ray* ray, Triangle* t) {
    Vector* e1 = new Vector(subtract(t->vertices[1], t->vertices[0]));
    Vector* e2 = new Vector(subtract(t->vertices[2], t->vertices[0]));
    Vector* cross = crossProduct(ray->direction, e2);
    float det = dotProduct(cross, e1);
    if ((det > -0.00001) && (det < 0.00001)) {
        return 0.0; //no intersection
    }
    Vector* toOrigin = new Vector(subtract(ray->origin, t->vertices[0]));
    float u = dotProduct(toOrigin, cross) / det;
    if ((u < 0.0) || (u > 1.0)) {
        return 0.0; //outside the triangle
    }
    Vector* cross1 = crossProduct(toOrigin, e1);
    float v = dotProduct(ray->direction, cross1) / det;
    if ((v < 0.0) || (v + u > 1.0)) {
        return 0.0; //outside the triangle
    }
    return dotProduct(e2, cross1) / det;
}

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
    float* zBuff;
    std::vector<Sphere*> spheres;
    std::vector<Triangle*> triangles;
    std::map<int, int> map;
    Scene(Point* ll1, Point* lr1, Point* ul1, Point* ur1, Point* eye1, int w, int h) {
        width = w;
        height = h;
        depth = -1*(pow(w*w + h*h, 0.5))*((ll1->z + lr1->z + ul1->z + ur1->z)/4);
        //printf("%d\n", depth);
        ll = new Point(ll1->x * w/2, ll1->y * h/2, ll1->z * depth);
        lr = new Point(lr1->x * w/2, lr1->y * h/2, lr1->z * depth);
        ul = new Point(ul1->x * w/2, ul1->y * h/2, ul1->z * depth);
        ur = new Point(ur1->x * w/2, ur1->y * h/2, ur1->z * depth);
        eye = new Point(eye1->x * w/2, eye1->y * h/2, eye1->z * depth);
        rsamples = new float[w*h];
        gsamples = new float[w*h];
        bsamples = new float[w*h];
        zBuff = new float[w*h];
    }
    void add_sphere(Sphere* s) {
        spheres.push_back(s);
        map[s->sid] = 1;
    }
    void add_tri(Triangle* t) {
        triangles.push_back(t);
        map[t->sid] = 0;
    }
    bool blocked(Ray* ray, int sid, Shape** blocker) {
        int i;
        float bt = FLT_MAX;
        for (i=0;i<spheres.size();i++) {
            //printf("spheres\n");
            Sphere* s = spheres.at(i);
            if (s->sid != sid) {
                float t = time_to_sphere(ray, s);
                //printf("%f %d %d\n", t, sid, s->sid);
                if (t > 0.0001) {
                    if (t < bt) {
                        //printf("small\n");
                        bt = t;
                        //printf("%f, %f \n", bt, FLT_MAX);
                        *blocker = s;
                    }
                }
            }
        }
        for (i=0;i<triangles.size();i++) {
            //printf("tri\n");
            Triangle* t1 = triangles.at(i);
            if (t1->sid != sid) {
                float t = time_to_triangle(ray, t1);
                //printf("%f\n", t);
                if (t > 0.0001) {
                    if (t < bt) {
                        bt = t;
                        *blocker = t1;
                    }
                }
            }
        }
        //printf("done\n");
        if ((*blocker) != NULL) {
            //printf("hit\n");
            //printf("%d\n", (*blocker)->sid);
            return true;
        }
        return false;
    }
};

Color* get_color(Ray* ray, Scene* my_scene, int index, int ttl) {
    float r = 0.0;
    float g = 0.0;
    float b = 0.0;
    if (ttl > 0) {
        Shape* s1 = NULL;
        //printf("test\n");
        if (my_scene->blocked(ray, 0, &s1)) {
            //printf("hit\n");
            if (my_scene->map[s1->sid] == 1) {
                Sphere* s = (Sphere*) s1;
                float t = time_to_sphere(ray, s);
                if (t>0) {
                    Point* worldIntersect = make_point(scale(ray->direction, t), ray->origin);
                    Point* intersect = normalize(scalez(add(worldIntersect, s->center), -1));
                    Vector* normal = normalize(new Vector(intersect));
                    Vector* view = new Vector(0,0,1);
                    int i;
                    for(i=0;i<plights.size();i++) {
                        Light* light = plights.at(i);
                        Vector* ltos = normalize(new Vector(scale(light->pos, s->radius), intersect));
                        Ray* shadowRay = new Ray(worldIntersect, normalize(new Vector(scale(light->pos, s->radius), worldIntersect)));
                        Shape* blocker = NULL;
                        if (!my_scene->blocked(shadowRay, s->sid, &blocker)) {
                            Vector* reflection = reflect(normalize(ltos), normal);
                            float diff = max(dotProduct(normalize(ltos), normal), zero);
                            float spec = pow(max(dotProduct(view, normalize(reflection)), zero), s->p);
                            Ray* reflect = new Ray(worldIntersect, normalize(reflection));
                            Shape* reflecter = NULL;
                            if (my_scene->blocked(reflect, s->sid, &reflecter)) {
                                Color* c = get_color(reflect, my_scene, index, ttl-1);
                                r = r + s->kr->r*c->r;
                                g = g + s->kr->g*c->g;
                                b = b + s->kr->b*c->b;
                            }
                            r = r + s->ka->r*light->color->r + s->kd->r*light->color->r*diff + s->ks->r*light->color->r*spec;
                            g = g + s->ka->g*light->color->g + s->kd->g*light->color->g*diff + s->ks->g*light->color->g*spec;
                            b = b + s->ka->b*light->color->b + s->kd->b*light->color->b*diff + s->ks->b*light->color->b*spec;
                        } else {
                            //printf("blocked\n");
                            r = r + s->ka->r*light->color->r;
                            g = g + s->ka->g*light->color->g;
                            b = b + s->ka->b*light->color->b;
                        }
                    }
                    for(i=0;i<dlights.size();i++) {
                        Light* light = dlights.at(i);
                        Vector* negLight = new Vector(scale(light->pos, -1));
                        Ray* shadowRay = new Ray(worldIntersect, normalize(new Vector(scale(light->pos, -1), worldIntersect)));
                        //printf("%f %f %f\n", negLight->x, negLight->y, negLight->z);
                        //printf("%f %f %f\n", shadowRay->direction->x, shadowRay->direction->y, shadowRay->direction->z);
                        Shape* blocker = NULL;
                        if (!my_scene->blocked(shadowRay, s->sid, &blocker)) {
                            Vector* reflection = reflect(normalize(negLight), normal);
                            float diff = max(dotProduct(normalize(negLight), normal), zero);
                            float spec = pow(max(dotProduct(view, normalize(reflection)), zero), s->p);
                            Ray* reflect = new Ray(worldIntersect, normalize(reflection));
                            Shape* reflecter = NULL;
                            if (my_scene->blocked(reflect, s->sid, &reflecter)) {
                                Color* c = get_color(reflect, my_scene, index, ttl-1);
                                r = r + s->kr->r*c->r;
                                g = g + s->kr->g*c->g;
                                b = b + s->kr->b*c->b;
                            }
                            r = r + s->ka->r*light->color->r + s->kd->r*light->color->r*diff + s->ks->r*light->color->r*spec;
                            g = g + s->ka->g*light->color->g + s->kd->g*light->color->g*diff + s->ks->g*light->color->g*spec;
                            b = b + s->ka->b*light->color->b + s->kd->b*light->color->b*diff + s->ks->b*light->color->b*spec;
                        } else {
                            //printf("blocked\n");
                            r = r + s->ka->r*light->color->r;
                            g = g + s->ka->g*light->color->g;
                            b = b + s->ka->b*light->color->b;
                        }

                    }
                }
            } else {
                Triangle* s = (Triangle*) s1;
                float t = time_to_triangle(ray, s);
                if (t>0.0) {
                    Point* worldIntersect = make_point(scale(ray->direction, t), ray->origin);
                    Point* intersect = scalez(add(worldIntersect, s->center), -1);
                    Vector* normal = normalize(new Vector(intersect));
                    Vector* view = new Vector(0,0,1);
                    int i;
                    for(i=0;i<plights.size();i++) {
                        Light* light = plights.at(i);
                        Vector* ltos = normalize(new Vector(scale(light->pos, s->radius), intersect));
                        Ray* shadowRay = new Ray(worldIntersect, normalize(new Vector(scale(light->pos, s->radius), worldIntersect)));
                        Shape* blocker = NULL;
                        if (!my_scene->blocked(shadowRay, s->sid, &blocker)) {
                            Vector* reflection = reflect(normalize(ltos), normal);
                            float diff = max(dotProduct(normalize(ltos), normal), zero);
                            float spec = pow(max(dotProduct(view, normalize(reflection)), zero), s->p);
                            Ray* reflect = new Ray(worldIntersect, normalize(reflection));
                            Shape* reflecter = NULL;
                            if (my_scene->blocked(reflect, s->sid, &reflecter)) {
                                Color* c = get_color(reflect, my_scene, index, ttl-1);
                                r = r + s->kr->r*c->r;
                                g = g + s->kr->g*c->g;
                                b = b + s->kr->b*c->b;
                            }
                            r = r + s->ka->r*light->color->r + s->kd->r*light->color->r*diff + s->ks->r*light->color->r*spec;
                            g = g + s->ka->g*light->color->g + s->kd->g*light->color->g*diff + s->ks->g*light->color->g*spec;
                            b = b + s->ka->b*light->color->b + s->kd->b*light->color->b*diff + s->ks->b*light->color->b*spec;
                        } else {
                            //printf("blocked\n");
                            r = r + s->ka->r*light->color->r;
                            g = g + s->ka->g*light->color->g;
                            b = b + s->ka->b*light->color->b;
                        }
                    }
                    for(i=0;i<dlights.size();i++) {
                        Light* light = dlights.at(i);
                        Vector* negLight = new Vector(scale(light->pos, -1 * s->radius));
                        Ray* shadowRay = new Ray(worldIntersect, normalize(new Vector(scale(light->pos, -1*s->radius), worldIntersect)));
                        Shape* blocker = NULL;
                        if (!my_scene->blocked(new Ray(worldIntersect, normalize(negLight)), s->sid, &blocker)) {
                            Vector* reflection = reflect(normalize(negLight), normal);
                            float diff = max(dotProduct(normalize(negLight), normal), zero);
                            float spec = pow(max(dotProduct(view, normalize(reflection)), zero), s->p);
                            Ray* reflect = new Ray(worldIntersect, normalize(reflection));
                            Shape* reflecter = NULL;
                            if (my_scene->blocked(reflect, s->sid, &reflecter)) {
                                Color* c = get_color(reflect, my_scene, index, ttl-1);
                                r = r + s->kr->r*c->r;
                                g = g + s->kr->g*c->g;
                                b = b + s->kr->b*c->b;
                            }
                            r = r + s->ka->r*light->color->r + s->kd->r*light->color->r*diff + s->ks->r*light->color->r*spec;
                            g = g + s->ka->g*light->color->g + s->kd->g*light->color->g*diff + s->ks->g*light->color->g*spec;
                            b = b + s->ka->b*light->color->b + s->kd->b*light->color->b*diff + s->ks->b*light->color->b*spec;
                        } else {
                            //printf("blocked\n");
                            r = r + s->ka->r*light->color->r;
                            g = g + s->ka->g*light->color->g;
                            b = b + s->ka->b*light->color->b;
                        }

                    }
                }
            }
        }
    }
    return new Color(r,g,b);
}


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
                Point* pos = new Point(i - width/2, (j - height/2)*-1, depth);
                //printf("About to make Vector\n");
                Ray* ray = new Ray(my_scene->eye, new Vector(pos, my_scene->eye));
                //printf("About to Generate\n");
                shoot_ray(ray, i + j*width);
            }
        }
    }
    void shoot_ray(Ray* ray, int index) {
        Color* c = get_color(ray, my_scene, index, 2);
        my_scene->rsamples[index] = c->r;
        my_scene->gsamples[index] = c->g;
        my_scene->bsamples[index] = c->b;
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
               fprintf(file, " %d %d %d ", r, g, b);
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
    int w = 500;
    int h = 500;
    Point* ll = new Point(-1,-1,-3);
    Point* lr = new Point(1,-1,-3);
    Point* ul = new Point(-1,1,-3);
    Point* ur = new Point(1,1,-3);
    Point* eye = new Point(0,0,0);
    //printf("Start\n");
    //printf("%d %d %d\n", argc, i, i<argc);
    Scene* scene = new Scene(ll, lr, ul, ur, eye, w, h);
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
        // else if (!(strcmp(argv[i],"-f"))) {
        //     h = fopen(argv[i+1], "rb");
        //     i += 2;
        // }
    }
    Color* ka = new Color(0.1, 0.1, 0.1);
    Sphere* s1 = new Sphere(0,0,-20,3, ka, new Color(1, 0, 1), new Color(1,1,1), 50, new Color(0,0,0), 1);
    Sphere* s2 = new Sphere(-2,2,-15,1, ka, new Color(1, 1, 0), new Color(1,1,1), 50, new Color(0,0,0), 2);
    //Sphere* s3 = new Sphere(0,-1,-2,0.25, ka, new Color(0, 1, 0), new Color(1,1,1), 50, new Color(1,1,1), 3);
    Sphere* s4 = new Sphere(-2,-2,-15,1, ka, new Color(0, 1, 1), new Color(1,1,1), 50, new Color(0,0,0), 4);
    //Sphere* s5 = new Sphere(1,0,-2,0.25, ka, new Color(1, 0, 1), new Color(1,1,1), 50, new Color(1,1,1), 5);
    Triangle* t = new Triangle(new Point(5, 5, -17), new Point(1, 4, -20), new Point(-6, -1, -20), ka, new Color(0.1, 0.1, 0.8), new Color(1, 1, 1), 50, new Color(1,1,1), 3);
    scene->add_sphere(s2);
    scene->add_sphere(s1);
    //scene->add_sphere(s3);
    scene->add_sphere(s4);
    //scene->add_sphere(s5);
    scene->add_tri(t);
    //printf("lights\n");
    //plights.push_back(new Light(new Point(-2, 2, 2), new Color(1, 1, 1)));
    //plights.push_back(new Light(new Point(-2, -2, 2), new Color(.6, .6, .6)));
    //plights.push_back(new Light(new Point(-2, 2, 2), new Color(.6, .6, .6)));
    dlights.push_back(new Light(new Point(2, -2, -2), new Color(1, 1, 1)));
    dlights.push_back(new Light(new Point(2, 2, -2), new Color(0, 0, 1)));
    //printf("sampler\n");
    Sampler* sampler = new Sampler(scene);
    Film* my_film = new Film(sampler);
    my_film->writeImage();
    return 0;
}