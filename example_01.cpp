
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

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>
#include <string.h>


#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }

using namespace std;

//****************************************************
// Some Classes
//****************************************************

class Viewport;

class Viewport {
  public:
    int w, h; // width and height
};

struct Color {
  float r, g, b;
  Color(float r1=0.0, float g1=0.0, float b1=0.0) {
    r = r1;
    g = g1;
    b = b1;
  }
}; 

struct Vector {
  float x, y, z;
  float mag;
  Vector* unit;
  Vector(float x1, float y1, float z1) {
    x = x1;
    y = y1;
    z = z1;
    mag = pow(sqr(x1) + sqr(y1) + sqr(z1), 0.5);
    if (mag != 1) {
      unit = new Vector(x1/mag, y1/mag, z1/mag);
    }
  }
};

struct Light {
  Vector* pos;
  Color* color;
  Light(float x1, float y1, float z1, float r, float g, float b) {
    pos = new Vector(x1, y1, z1);
    color = new Color(r, g, b);
  }
};

//****************************************************
// Global Variables
//****************************************************
Viewport	viewport;

int p = 0; //Phong
int pointLights = 0;
Light* pLights[5];
int directionalLights = 0;
Light* dLights[5];
Color* ka;
Color* kd;
Color* ks;
Vector* viewer = new Vector(0,0,1);
int id;


//****************************************************
// Simple init function
//****************************************************
void initScene(){

  // Nothing to do here for this simple example.

}


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glViewport (0,0,viewport.w,viewport.h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, viewport.w, 0, viewport.h);

}


//****************************************************
// A routine to set a pixel by drawing a GL point.  This is not a
// general purpose routine as it assumes a lot of stuff specific to
// this example.
//****************************************************

void setPixel(int x, int y, GLfloat r, GLfloat g, GLfloat b) {
  glColor3f(r, g, b);
  glVertex2f(x + 0.5, y + 0.5);   // The 0.5 is to target pixel
  // centers 
  // Note: Need to check for gap
  // bug on inst machines.
}

float dotProduct(Vector a, Vector b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector* subtract(Vector a, Vector b) {
  return new Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector* add(Vector a, Vector b) {
  return new Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector* scale(Vector a, float s) {
  return new Vector(a.x * s, a.y * s, a.z * s);
}

Vector* reflect(Vector l, Vector n) {
  float dot = dotProduct(l, n);
  return add(*scale(*scale(n, dot), 2), *scale(l, -1));
}

//****************************************************
// Draw a filled circle.  
//****************************************************


void circle(float centerX, float centerY, float radius) {
  // Draw inner circle
  glBegin(GL_POINTS);

  // We could eliminate wasted work by only looping over the pixels
  // inside the sphere's radius.  But the example is more clear this
  // way.  In general drawing an object by loopig over the whole
  // screen is wasteful.

  int i,j;  // Pixel indices

  int minI = max(0,(int)floor(centerX-radius));
  int maxI = min(viewport.w-1,(int)ceil(centerX+radius));

  int minJ = max(0,(int)floor(centerY-radius));
  int maxJ = min(viewport.h-1,(int)ceil(centerY+radius));



  for (i=0;i<viewport.w;i++) {
    for (j=0;j<viewport.h;j++) {

      // Location of the center of pixel relative to center of sphere
      float x = (i+0.5-centerX);
      float y = (j+0.5-centerY);

      float dist = sqrt(sqr(x) + sqr(y));

      if (dist<=radius) {
        float r = 0.0;
        float g = 0.0;
        float b = 0.0;

        // This is the front-facing Z coordinate
        float z = sqrt(radius*radius-dist*dist);
        Vector* surface = new Vector(x, y, z);

        for (int pl=0;pl<pointLights;pl++) {
          Light* light = pLights[pl];
          Vector* ltos = subtract(*light->pos, *surface);
          float intensity = ltos->mag;
          Vector* reflection = reflect(*ltos, *surface);
          r = r + ka->r * intensity + kd->r * intensity * dotProduct(*ltos->unit, *surface->unit)
              + ks->r * intensity * pow(dotProduct(*viewer, *reflection->unit), p);
          g = g + ka->g * intensity + kd->g * intensity * dotProduct(*ltos->unit, *surface->unit)
              + ks->g * intensity * pow(dotProduct(*viewer, *reflection->unit), p);
          b = b + ka->b * intensity + kd->b * intensity * dotProduct(*ltos->unit, *surface->unit)
              + ks->b * intensity * pow(dotProduct(*viewer, *reflection->unit), p);
        }
        for (int dl=0;dl<directionalLights;dl++) {
          Light* light = dLights[dl];
          Vector* negLight = scale(*light->pos, -1);
          float intensity = light->pos->mag;
          Vector* reflection = reflect(*negLight, *surface);
          r = r + ka->r * intensity + kd->r * intensity * dotProduct(*negLight->unit, *surface->unit)
              + ks->r * intensity * pow(dotProduct(*viewer, *reflection->unit), p);
          g = g + ka->g * intensity + kd->g * intensity * dotProduct(*negLight->unit, *surface->unit)
              + ks->g * intensity * pow(dotProduct(*viewer, *reflection->unit), p);
          b = b + ka->b * intensity + kd->b * intensity * dotProduct(*negLight->unit, *surface->unit)
              + ks->b * intensity * pow(dotProduct(*viewer, *reflection->unit), p);
        }
        setPixel(i,j, r, g, b);

        // This is amusing, but it assumes negative color values are treated reasonably.
        // setPixel(i,j, x/radius, y/radius, z/radius );
      }


    }
  }


  glEnd();
}
//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {

  glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer

  glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
  glLoadIdentity();				        // make sure transformation is "zero'd"


  // Start drawing
  circle(viewport.w / 2.0 , viewport.h / 2.0 , min(viewport.w, viewport.h) / 3.0);

  glFlush();
  glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

void myKey(unsigned char key, int x, int y) {
  if (key == 32) {
    glutDestroyWindow(id);
  }
}

//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {
  //This initializes glut
  glutInit(&argc, argv);

  //This tells glut to use a double-buffered window with red, green, and blue channels 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  // Initalize theviewport size
  viewport.w = 400;
  viewport.h = 400;

  //The size and position of the window
  glutInitWindowSize(viewport.w, viewport.h);
  glutInitWindowPosition(0,0);
  id = glutCreateWindow(argv[0]);

  int i = 1;
  while (i < argc) {
    if (!(strcmp(argv[i],"-ka"))) {
      ka = new Color(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
      i += 4;
    }
    if (!(strcmp(argv[i],"-ks"))) {
      ks = new Color(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
      i += 4;
    }
    if (!(strcmp(argv[i],"-kd"))) {
      kd = new Color(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
      i += 4;
    }
    if (!(strcmp(argv[i],"-pl"))) {
      if (pointLights < 5) {
        pLights[pointLights] = new Light(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]),
                                          atof(argv[i+4]), atof(argv[i+5]), atof(argv[i+6]));
        pointLights += 1;
        i += 7;
      }
    }
    if (!(strcmp(argv[i],"-dl"))) {
      if (directionalLights < 5) {
        dLights[directionalLights] = new Light(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]),
                                          atof(argv[i+4]), atof(argv[i+5]), atof(argv[i+6]));
        directionalLights += 1;
        i += 7;
      }
    }
    if (!(strcmp(argv[i],"-sp"))){
      p = atof(argv[i+1]);
      i += 2;
    }
  }

  initScene();							// quick function to set up scene

  glutDisplayFunc(myDisplay);				// function to run when its time to draw something
  glutReshapeFunc(myReshape);				// function to run when the window gets resized
  glutKeyboardFunc(myKey);

  glutMainLoop();							// infinite loop that will keep drawing and resizing
  // and whatever else

  return 0;
}








