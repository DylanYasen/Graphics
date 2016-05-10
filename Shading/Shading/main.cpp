//
//  main.cpp
//  Shading
//
//  Created by Yadikaer Yasheng on 5/2/16.
//  Copyright © 2016 Yadikaer Yasheng. All rights reserved.
//

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include "MVectorMatrix.h"

using namespace __svl_lib;

#define EPSILON 0.000001f

struct Object {
    std::string name;
    Vector3f center;
    Vector3f color;
    Vector3f normal;
    Vector3f diffuseColor;
    Vector3f specularColor;
    virtual ~Object(){}
    virtual float RayIntersection(Vector3f, Vector3f)=0;
    virtual Vector3f GetNormal(Vector3f){return normalize(normal);}
};

struct Sphere : public Object {
    float radius;
    
    float RayIntersection(Vector3f rayDir, Vector3f rayOrigin) override{
        Vector3f d = rayDir;
        Vector3f p = rayOrigin;
        
        float A = d.dot(d);
        float B = 2 * d.dot(p - center);
        float C = (p - center).dot(p - center) - radius * radius;
        float delta = B * B - 4 * A * C;
        
        if(delta < 0){
            return INFINITY;
        }
        float del = sqrtf(delta);
        float t1 = (-B + del) / (2 * A);
        float t2 = (-B - del) / (2 * A);
        
        if(t1 < 0 && t2 < 0)
            return INFINITY;
        else if (t1 < 0 && t2 > 0)
            return t2;
        else if (t2 < 0 && t1 > 0 )
            return t1;
        else if (t1 >= 0 && t2 >= 0)
            return fminf(t1, t2);
        
        return INFINITY;
    }
    
    Vector3f GetNormal(Vector3f point) override{
        return normalize(center - point);
    }
};

struct Plane : public Object{
    
    float RayIntersection(Vector3f rayDir, Vector3f rayOrigin) override{
    
        Vector3f p0 = center;
        Vector3f d = normalize(rayDir);
        Vector3f e = rayOrigin;
        Vector3f n = GetNormal(Vector3f(0,0,0));
        
        /*
         Ray: P = e + td
         Plane : (P−P0)⋅n = 0
         Subsitute for P:
         (e + td − P0)⋅n = 0
         (e - P0)⋅n = -td⋅n
         t = ((e - P0)⋅n) / (-d⋅n)
        */
        
        float denom = (d).dot(n);
        
        // ray parallel to the plane
        if (denom == 0.0) {
            printf("failed\n");
            return INFINITY;
        }
        
        float num = (p0 - e).dot(n);
        float t = num / denom;
        
        if(t < 0)
            return INFINITY;
        
        return t;
    }
};

struct Scene {
    std::string name;
    std::vector<Object*> objects;
    
    // view
    Vector3f viewPoint;
    Vector3f viewDir;
    Vector3f viewUp;
    float projDistance;
    float viewWidth;
    float viewHeight;
    float imageSize;
    
    // lighting
    float specularExponent;
    Vector3f lightPosition;
    Vector3f lightIintensity;

    
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

Vector3f viewPoint(0,0,5);
Vector3f viewDir(0,0,-5);
Vector3f viewUp(0,1,0);
float projDistance = 5;
float viewWidth = 2.5;
float viewHeight = 2.5;
float imageSize = 300;

float specularExponent = 80;
Vector3f lightPosition(6,0,10);
Vector3f lightIintensity(0.2,0.2,0.2);
float Ka = 0.5;
float Kd  = 0.5;
float Ks = 0.5;

Vector3f bgColor(1,1,1);

bool CheckShadow(Scene scene,Vector3f hitPoint, Vector3f lightRayDir){
    
    // avoid numeric imprecision
    hitPoint += normalize(lightRayDir) * 0.2;
    
    // check if lightray is blocked by scene objects
    for (int i = 0; i < scene.objects.size(); i++) {
        float t = scene.objects[i]->RayIntersection(lightRayDir,hitPoint);
        if (t != INFINITY) {
            //printf("%f\n",t);
            printf("shadow ray hit %s \n",scene.objects[i]->name.c_str());
            return true;
        }
    }
    return false;
}

Vector3f Shading(Object *obj, Vector3f hitPoint, Vector3f viewRay){
  
    // v+l h= ∥v+l∥,
    // L = kd I max(0, n · l) + ks I max(0, n · h)^p,
    Vector3f n = obj->GetNormal(hitPoint);
    Vector3f l = normalize(hitPoint - lightPosition);
    Vector3f v = -normalize(viewRay);
    Vector3f h = normalize(v + l);
    
    Vector3f diffuse = Kd * obj->diffuseColor * fmaxf(0, n.dot(l));
    Vector3f specular = Ks * obj->specularColor * pow(fmaxf(0, n.dot(h)), specularExponent);
    Vector3f ambient = Ka * lightIintensity;
 
    // light ray not blocked
    Vector3f L = diffuse + specular + ambient;
    
    // clamp 0 - 1
    if(L[0] > 1)
        L[0] = 1;
    if(L[1] > 1)
        L[1] = 1;
    if(L[2] > 1)
        L[2] = 1;
   
    
    //obj->color = L;
    return L;
}

void WriteToPPM(Vector3f imageBuffer[][300],int w, int h, std::string fn){
    
    FILE* fp;
    if (NULL != (fp = fopen(fn.c_str(), "w"))){
        // Write the 'header' information
        fprintf(fp, "P3\n");
        fprintf(fp, "%d %d %d\n",w,h,255);
        
        for (int r = h-1; r >= 0; r--) {
            for (int c = 0; c < w; c++) {
                
                // scale to 255
                imageBuffer[r][c] *= 255;
                
                //printf("%f %f %f\n",imageBuffer[r][c][0],imageBuffer[r][c][1],imageBuffer[r][c][2]);
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
    
    printf("%lu\n",scene.objects.size());
    
    // establish camera basis
    Vector3f w = -normalize(viewDir);
    Vector3f u = normalize(w.cross(viewUp));
    Vector3f v = normalize(u.cross(w));
    
//    printf("u: %f %f %f\n",u[0],u[1],u[2]);
//    printf("v: %f %f %f\n",v[0],v[1],v[2]);
//    
    
    // unit converion
    float widthPerPx = viewWidth / imageSize;
    float heightPerPx = viewHeight / imageSize;
    
    Vector3f e = viewPoint;
    
    // find bottom left point & view dir
    Vector3f eyeBotLeft = e - (viewWidth / 2) * u - (viewHeight / 2) * v + (widthPerPx / 2) * u + (heightPerPx / 2) * v;
    
//    printf("bot left :%f %f %f\n",eyeBotLeft[0],eyeBotLeft[1],eyeBotLeft[2]);
    
    // shoot ray from each pixel position
    for (int r = 0; r < imageSize; r++) {
        for (int c = 0; c < imageSize; c++) {
            
            Vector3f s = eyeBotLeft + u * c * widthPerPx + v * r * heightPerPx;
            Vector3f p = s;
            Vector3f d = -w ;
            
            //printf("%f %f %f\n",p[0],p[1],p[2]);
            
            for (int i = 0; i < scene.objects.size(); i++) {
                
                float t = scene.objects[i]->RayIntersection(d, p);
                
                // ray hits object
                if (t != INFINITY) {
                    Vector3 hitpoint = p + d * t;
                    Vector3 viewRay = p - hitpoint;
                    
                    // TODO: move ambient to scene object
                    Vector3f ambient = Ka * lightIintensity;
                    Vector3f lightray = lightPosition - hitpoint;
                    
                    if(CheckShadow(scene, hitpoint, lightray)){
                        imageBuffer[r][c] = ambient;
                    }
                    else{
                        Vector3 color = Shading(scene.objects[i], hitpoint, viewRay);
                        imageBuffer[r][c] = color;
                    }
                }
                else{
//                    printf("ray doesn't hit %s \n",scene.objects[i]->name.c_str());
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
    scene1.name = "0-three-sphere";
    
    Plane *plane = new Plane();
    plane->center = Vector3f(0,0,0);
    plane->normal = Vector3f(0,0,-1);
    plane->diffuseColor = Vector3f(0,1,1);
    plane->specularColor = Vector3f(1,1,1);
    plane->name = "plane";
    scene1.objects.push_back(plane);

    Sphere *sphere0 = new Sphere();
    sphere0->center = Vector3f(-0.5,0,0.4);
    sphere0->radius = 0.4f;
    sphere0->diffuseColor = Vector3f(1,1,0);
    sphere0->specularColor = Vector3f(1,1,1);
    sphere0->name = "sphere0";
 
    Sphere *sphere1 = new Sphere();
    sphere1->center = Vector3f(0,0.4,0.4);
    sphere1->radius = 0.4f;
    sphere1->diffuseColor = Vector3f(0.5,1,0.2);
    sphere1->specularColor = Vector3f(1,1,1);
    sphere1->name = "sphere1";
   
    Sphere *sphere2 = new Sphere();
    sphere2->center = Vector3f(0.5,0,0.4);
    sphere2->radius = 0.4f;
    sphere2->diffuseColor = Vector3f(0.6,0.2,0.7);
    sphere2->specularColor = Vector3f(1,1,1);
    sphere2->name = "sphere2";
    
    // kind of a render queue
    scene1.objects.push_back(sphere1);
    scene1.objects.push_back(sphere0);
    scene1.objects.push_back(sphere2);
    
    ResetImageBuffer(imageBuffer);
    Trace(scene1,imageBuffer);
    std::string fp = "/Users/Yadikaer/Projects/Graphics/Shading/Shading/"+ scene1.name+ "image.ppm";
    WriteToPPM(imageBuffer, imageSize, imageSize, fp);
    scene1.Clear();
    //
    //  ====================================================
    //
    
    //
    //  ================== Scene 2 ==================
    //
    /*
    Scene scene2;
    scene2.name = "1-multiple-spheres";
    sphere0 = new Sphere();
    sphere0->center = Vector3f(0,0,0);
    sphere0->radius = 0.4;
    sphere0->color = Vector3f(1,0,0);
    scene2.objects.push_back(sphere0);
    
    sphere1 = new Sphere();
    sphere1->center = Vector3f(0,1,0);
    sphere1->radius = 0.4;
    sphere1->color = Vector3f(0,1,0);
    scene2.objects.push_back(sphere1);
    
    sphere2 = new Sphere();
    sphere2->center = Vector3f(0,-1,0);
    sphere2->radius = 0.4;
    sphere2->color = Vector3f(0,0,1);
    scene2.objects.push_back(sphere2);
    
    ResetImageBuffer(imageBuffer);
    Trace(scene2,imageBuffer);
    fp = "/Users/Yadikaer/Projects/Graphics/Shading/Shading/"+ scene2.name+ "image.ppm";
    WriteToPPM(imageBuffer, imageSize, imageSize, fp);
    scene2.Clear();
     */
    
    //  ====================================================

    return 0;
}
