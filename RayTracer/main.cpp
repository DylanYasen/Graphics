//
//  main.cpp
//  RayTracer
//
//  Created by Yadikaer Yasheng on 4/24/16.
//  Copyright © 2016 Yadikaer Yasheng. All rights reserved.
//

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include "MVectorMatrix.h"

using namespace __svl_lib;

#define EPSILON 0.000001

struct Object {
    Vector3f center;
    Vector3f color;
    virtual ~Object(){}
    virtual float RayIntersection(Vector3f, Vector3f)=0;
};

struct Sphere : public Object {
    float radius;
   
    float RayIntersection(Vector3f d, Vector3f p) override{
        float A = d.dot(d);
        float B = 2 * d.dot(p - center);
        float C = (p - center).dot(p - center) - radius * radius;
        float delta = B * B - 4 * A * C;
        
        float d = sqrtf(b*b)
        float t1 = -b + sqrtf(b*b - )
        
        return delta;
    }
};

struct Triangle : public Object {
    
    Vector3f va,vb,vc;
    
    float RayIntersection(Vector3f d, Vector3f p) override{
        
        // calculate vertices in world coordinates
        Vector3f A,B,C;
        A = center + va;
        B = center + vb;
        C = center + vc;
        
        // edges
        Vector3f AB = B - A;
        Vector3f AC = C - A;
        Vector3f BC = C - B;
        
        // find tri plane normal
        Vector3f n = AC.cross(AB);
        n.normalize();
        
        // check if on ray
        float t = (A.dot(n) - p.dot(n))/d.dot(n);
        
        // parallel to the plane
        if (t == 0) {
            return -1;
        }
       
        // point of intersection on the plane
        Vector3f I = p + d*t;
        
        float alpha = AB.cross((I-A)).dot(n);
        float beta = BC.cross((I-B)).dot(n);
        float gamma = AC.cross((I-C)).dot(n);
                               
        if(alpha > 0 && beta > 0 && gamma > 0)
            return 1;
        
        return -1;
        
        // check if inside of the triangle
    }
};

struct Slab : public Object {
    float xmin,ymin,xmax,ymax;
    
    float RayIntersection(Vector3f d, Vector3f p) override{
    
        float txmin = (xmin - p[0]) / d[0];
        float txmax = (xmax - p[0]) / d[0];
        float tymin = (ymin - p[1]) / d[1];
        float tymax = (ymax - p[1]) / d[1];
        
        float tmin = std::max(txmin,tymin);
        float tmax = std::min(txmax,tymax);
        
        if(tmin > tymax || tymin > tmax)
            return -1;
        
        return 1;
    }
    
};

struct Box : public Object {
    Vector3f size;
    
    float RayIntersection(Vector3f d, Vector3f p) override{
        
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        
        if (d[0] >= 0) {
            tmin = ((center[0] - size[0]) - p[0]) / d[0];
            tmax = ((center[0] + size[0]) - p[0]) / d[0];
        }
        else{
            tmin = ((center[0] + size[0]) - p[0]) / d[0];
            tmax = ((center[0] - size[0]) - p[0]) / d[0];
        }
        
        if (d[1] >= 0) {
            tymin = ((center[1] - size[1])) / d[1];
            tymax = ((center[1] + size[1])) / d[1];
        }
        else{
            tymin = ((center[1] + size[1])) / d[1];
            tymax = ((center[1] - size[1])) / d[1];
        }
        
        if((tmin > tymax) || (tymin > tmax))
            return -1;
        
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;
        
        if (d[2] >= 0) {
            tzmin = (center[2] - d[2]) / d[2];
            tzmax = (center[2] + d[2]) / d[2];
        }
        else {
            tzmin = (center[2] + d[2]) / d[2];
            tzmax = (center[2] - d[2]) / d[2];
        }
        
        if ( (tmin > tzmax) || (tzmin > tmax) )
            return false;
        
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;
        
        return (1);
    }
};


struct Scene {
    std::string name;
    std::vector<Object*> objects;
    
    void Clear(){
        for(int i = 0; i < objects.size(); i++){
            delete (objects[i]);
        }
        objects.clear();
    }
};

void WriteToPPM(Vector3f imageBuffer[][300],int w, int h, std::string fn);
void Trace(Scene scene, Vector3f imageBuffer[][300]);
float RaySphereIntersection(Vector3f d, Vector3f p, Sphere sphere);
void ResetImageBuffer(Vector3f imageBuffer[][300]);


Vector3f viewPoint(5,4,3);
Vector3f viewDir(-5,-4,-3);
Vector3f projNormal(5,4,3);
Vector3f viewUp(0,1,0);
float projDistance = 5;
float viewWidth = 2.5;
float viewHeight = 2.5;
float imageSize = 300;

Vector3f diffuseColor (0.2, 0.3, 0.8);
Vector3f specularColor (1,1,0);
float specularExponent = 50;
Vector3f lightPosition(3,4,5);
Vector3f lightIintensity(1,1,1);
float Ka = 0.5;
float Kd  = 0.5;
float Ks = 0.5;

Vector3f bgColor(1,1,1);

Vector3f PhongShading(Vector3f hitPoint, Vector3f ray,Vector3f surfaceNormal){
    // v+l h= ∥v+l∥,
    // L = kd I max(0, n · l) + ks I max(0, n · h)^p,
    Vector3f n = surfaceNormal;
    Vector3f l = hitPoint - lightPosition;
    Vector3f v = ray;
    Vector3f h = normalize(v + l);
    
    Vector3f L = Kd * lightIintensity * fmaxf(0, n.dot(l)) + Ks * lightIintensity * pow(fmax(0, n.dot(h)), specularExponent);
    
    return L;
}

void WriteToPPM(Vector3f imageBuffer[][300],int w, int h, std::string fn){
    
    FILE* fp;
    if (NULL != (fp = fopen(fn.c_str(), "w"))){
        // Write the 'header' information
        fprintf(fp, "P3\n");
        fprintf(fp, "%d %d %d\n",w,h,1);
        
        for (int r = h-1; r >= 0; r--) {
            for (int c = 0; c < w; c++) {
                
                fprintf(fp, "%d %d %d\t",(int)imageBuffer[r][c][0],(int)imageBuffer[r][c][1],(int)imageBuffer[r][c][2]);
            }
            
            fprintf(fp, "\n");
        }
        
        fclose(fp);
    }
}

void ResetImageBuffer(Vector3f imageBuffer[][300]){
    
    for (int r = 0; r < 300; r++) {
        for (int c = 0; c < 300; c++) {
            imageBuffer[r][c] = Vector3(1,1,1);  // background color white
        }
    }
}

void Trace(Scene scene,Vector3f imageBuffer[][300]){
    
    // establish camera basis
    Vector3f w = -normalize(viewDir);
    Vector3f u = normalize(viewUp.cross(w));
    Vector3f v = normalize(w.cross(u));
    
    // unit converion
    float widthPerPx = viewWidth / imageSize;
    float heightPerPx = viewHeight / imageSize;

    Vector3f e = viewPoint;
    
    // find bottom left point & view dir
    Vector3f eyeBotLeft = e - (viewWidth / 2) * u - (viewHeight / 2) * v + Vector3f(widthPerPx/2,heightPerPx/2,0);
   
    // shoot ray from each pixel position
    for (int r = 0; r < imageSize; r++) {
        for (int c = 0; c < imageSize; c++) {
            
            Vector3f s = eyeBotLeft + u * c * widthPerPx + v * r * heightPerPx;
            Vector3f p = s;
            Vector3f d = -w ;
            
            for (int i = 0; i < scene.objects.size(); i++) {
                
                float delta = scene.objects[i]->RayIntersection(d, p);
                
                    if (delta >= 0) {
                        imageBuffer[r][c] = scene.objects[i]->color;
                    }
                }
            }
    }
}

int main(int argc, const char * argv[]) {
    
    Vector3f imageBuffer[300][300];
    
    //
    //  ================== Scene 1 ==================
    //
    Scene scene1;
    scene1.name = "0-single-unit-sphere";
    Sphere *sphere0 = new Sphere();
    sphere0->center = Vector3f(0,0,0);
    sphere0->radius = 1;
    sphere0->color = Vector3f(1,0,0);
    scene1.objects.push_back(sphere0);
    
    ResetImageBuffer(imageBuffer);
    Trace(scene1,imageBuffer);
    std::string fp = "/Users/Yadikaer/Projects/OpenGL/RayTracer/RayTracer/"+ scene1.name+ "image.ppm";
    WriteToPPM(imageBuffer, imageSize, imageSize, fp);
    scene1.Clear();
    //
    //  ====================================================
    //
    
    //
    //  ================== Scene 2 ==================
    //
    Scene scene2;
    scene2.name = "1-multiple-spheres";
    sphere0 = new Sphere();
    sphere0->center = Vector3f(0,0,0);
    sphere0->radius = 0.4;
    sphere0->color = Vector3f(1,0,0);
    scene2.objects.push_back(sphere0);

    Sphere *sphere1 = new Sphere();
    sphere1->center = Vector3f(0,1,0);
    sphere1->radius = 0.4;
    sphere1->color = Vector3f(0,1,0);
    scene2.objects.push_back(sphere1);
    
    Sphere *sphere2 = new Sphere();
    sphere2->center = Vector3f(0,-1,0);
    sphere2->radius = 0.4;
    sphere2->color = Vector3f(0,0,1);
    scene2.objects.push_back(sphere2);
    
    ResetImageBuffer(imageBuffer);
    Trace(scene2,imageBuffer);
    fp = "/Users/Yadikaer/Projects/OpenGL/RayTracer/RayTracer/"+ scene2.name+ "image.ppm";
    WriteToPPM(imageBuffer, imageSize, imageSize, fp);
    scene2.Clear();
    
    //  ====================================================
    
    
    //
    //  ================== Scene 3 ==================
    //
    Scene scene3;
    scene3.name = "2-triangle";
    Triangle *tri = new Triangle();
    tri->center = Vector3f(0,0,0);
    tri->va = Vector3f(-0.5,-0.5,0);
    tri->vb = Vector3f(0,0.5,0);
    tri->vc = Vector3f(0.5,-0.5,0);
    tri->color = Vector3f(0,0,1);
    scene3.objects.push_back(tri);
    
    ResetImageBuffer(imageBuffer);
    Trace(scene3,imageBuffer);
    fp = "/Users/Yadikaer/Projects/OpenGL/RayTracer/RayTracer/"+ scene3.name+ "image.ppm";
    WriteToPPM(imageBuffer, imageSize, imageSize, fp);
    scene3.Clear();
    
    //  ====================================================
    

    //
    //  ================== Scene 4 ==================
    //
    Scene scene4;
    scene4.name = "3-box";
    Box *box = new Box();
    box->size = Vector3(0.5,0.5,0.5);
    box->center = Vector3(0,0,0);
    box->color = Vector3(0,1,0);
    scene4.objects.push_back(box);
    
    ResetImageBuffer(imageBuffer);
    Trace(scene4,imageBuffer);
    fp = "/Users/Yadikaer/Projects/OpenGL/RayTracer/RayTracer/"+ scene4.name+ "image.ppm";
    WriteToPPM(imageBuffer, imageSize, imageSize, fp);
    scene4.Clear();
    
    //  ====================================================

    return 0;
}
