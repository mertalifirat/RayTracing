#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "cmath"
#include "pthread.h"
#define Epsilon 0.001
typedef unsigned char RGB[3];
using namespace std;
using namespace parser;


Scene scene;
unsigned char* image;
int w, h, lightNum;

typedef struct {
    int x, y, z;        
} vec3i;

typedef struct{
    Vec3f o, d;
    int height, width;
} Ray;

typedef struct{
    Camera currentCam;
    int beginI,beginJ;
    int endI,endJ;
    int h,w;

}ThreadStruct;

typedef enum {SPHERE, TRIANGLE, MESH} object_type;

struct Intersection {
    float hit_time;
    int material_id;
    vector<Face> faces;
    Face indices;
    int center_vertex_id;
    float radius;
    object_type type;
    Vec3f normal_vector;
    int index;
};


Vec3f scalarMult(const Vec3f& vec, float a){
    Vec3f result;
    result.x = vec.x * a;
    result.y = vec.y * a;
    result.z = vec.z * a;
    return result;
}


float dotProduct(const Vec3f& vec1,const Vec3f& vec2){
    return vec1.x*vec2.x+vec1.y*vec2.y+vec1.z*vec2.z;
}

Vec3f add(const Vec3f& vec1,const Vec3f& vec2){
    Vec3f result;
    result.x = vec1.x + vec2.x;
    result.y = vec1.y + vec2.y;
    result.z = vec1.z + vec2.z;
    return result;
}

Vec3f subtract(const Vec3f& vec1,const Vec3f& vec2){
    Vec3f result;
    result.x = vec1.x - vec2.x;
    result.y = vec1.y - vec2.y;
    result.z = vec1.z - vec2.z;
    return result;
}

Vec3f crossProduct(const Vec3f& vec1,const Vec3f& vec2){
    Vec3f result;
    result.x = vec1.y * vec2.z - vec1.z * vec2.y;
    result.y = vec1.z * vec2.x - vec1.x * vec2.z;
    result.z = vec1.x * vec2.y - vec1.y * vec2.x;
    return result;
}

Vec3f normalizeVector(const Vec3f& vec){
    Vec3f result;
    float length = sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z); 
    result.x = vec.x / length;
    result.y = vec.y / length;
    result.z = vec.z / length;
    return result;    
}

Ray generateRay(int i, int j, const Camera& camera){
    Ray ray;
    ray.o = camera.position;
    Vec3f gaze = camera.gaze, m, q;
    Vec3f unit_w = scalarMult(normalizeVector(camera.gaze), -1);
    Vec3f unit_v = normalizeVector(camera.up);
    Vec3f unit_u = normalizeVector(crossProduct(unit_v,unit_w));
    float l = camera.near_plane.x;
    float r = camera.near_plane.y;
    float b = camera.near_plane.z;
    float t = camera.near_plane.w;
    float distance = camera.near_distance;
    int w = camera.image_width, h = camera.image_height;
    m = add(ray.o, scalarMult(gaze, distance));
    q = add(m, add(scalarMult(unit_u, l), scalarMult(unit_v, t)));
    float s_u = (j + 0.5)*(r-l)/w; 
    float s_v = (i + 0.5)*(t-b)/h;
    Vec3f s = add(q, add(scalarMult(unit_u, s_u), scalarMult( unit_v, -1*s_v)));
    ray.d = normalizeVector(subtract(s, ray.o));
    ray.height = i;
    ray.width = j;
    return ray;
}


Intersection firstIntersection(const Intersection& hit1, const Intersection& hit2){
    return hit1.hit_time < hit2.hit_time ? hit1 : hit2 ;
}


Intersection sphere_hit(const Scene& scene, Ray ray){
    Intersection result;
    result.type = SPHERE;
    result.hit_time = MAXFLOAT;
    int min_sphere_index=0;
    float timeDiscr,r, time1, time2;
    Vec3f c,o=ray.o,d=ray.d;
    const vector<Sphere> &spheres = scene.spheres;
    const vector<Vec3f> &vertices = scene.vertex_data;
    bool hit_flag = false;
    for(int i = 0; i < spheres.size(); i++ ){
        Sphere sphere = spheres[i];
        r = sphere.radius;
        c = vertices[sphere.center_vertex_id - 1];
        Vec3f CE = subtract(ray.o, c);
        const float C = dotProduct(CE, CE) - r * r;
        const float B = 2 * dotProduct(ray.d, CE);
	    const float A = dotProduct(ray.d, ray.d);
	    const float timeDiscr = B*B - 4*A*C;
        if( timeDiscr >= 0){
            const float time1 = (-1 * B + sqrtf(timeDiscr))/2*A;
		    const float time2 = (-1 * B - sqrtf(timeDiscr))/2*A;
            if(fmin(time1, time2) >= 0 && result.hit_time > fmin(time1, time2)){result.hit_time = fmin(time1, time2);min_sphere_index = i;hit_flag = true;} 
        }
    }
    if(hit_flag){
        Sphere sphere = spheres[min_sphere_index];
        result.material_id = sphere.material_id;
        result.center_vertex_id = sphere.center_vertex_id;
        result.radius = sphere.radius;
        Vec3f intersection_point = add(ray.o,scalarMult(ray.d, result.hit_time));
        result.index = min_sphere_index;
        result.normal_vector = normalizeVector(scalarMult(subtract(intersection_point, vertices[sphere.center_vertex_id - 1]), (1/sphere.radius)));
    }
    return result;
} 


float magnitudeCalc(const Vec3f& vec){
    return sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);

}



Intersection triangle_hit(const Scene& scene, const Ray& ray){
    Intersection result;
    result.type = TRIANGLE;
    result.hit_time = MAXFLOAT;
    float t = MAXFLOAT; 
    const vector<Vec3f> &vertexData=scene.vertex_data;
    const vector<Triangle> &triangles=scene.triangles;
    const Vec3f &o=ray.o;
    const Vec3f &d=ray.d;
    int min_index = -1;
    bool hit_flag = false;
    for (int i=0; i<triangles.size(); i++){
        float tempT;
        int material_id=triangles[i].material_id;
        const Face& indices=triangles[i].indices;
        int v0=indices.v0_id,v1=indices.v1_id,v2=indices.v2_id;
        
        Vec3f a=vertexData[v0-1],b=vertexData[v1-1], c=vertexData[v2-1];
        tempT=dotProduct(subtract(a,o),indices.normal)/dotProduct(d,indices.normal);
        if (tempT>t){
            continue;
        }
        Vec3f hitPoint;
        hitPoint.x=o.x+d.x*tempT;
        hitPoint.y=o.y+d.y*tempT;
        hitPoint.z=o.z+d.z*tempT;
        Vec3f AC=subtract(c,a);
        float inv_AreaABC=1.0 / (magnitudeCalc(crossProduct(AC,indices.edge_2)));
        Vec3f PB=subtract(b,hitPoint),PC=subtract(c,hitPoint);
        float alfa=magnitudeCalc(crossProduct(PB,PC))*inv_AreaABC;
        
        Vec3f PA=subtract(a,hitPoint);
        float beta=magnitudeCalc(crossProduct(PC,PA))*inv_AreaABC;

        float sigma=magnitudeCalc(crossProduct(PA,PB))*inv_AreaABC;
        if(!(-Epsilon<=alfa && alfa<=1+Epsilon)){
            continue;
        }
        if(!(-Epsilon<=beta && beta<=1+Epsilon)){
            continue;
        }
        if(!(-Epsilon<=sigma && sigma<=1+Epsilon)){
            continue;
        }
        if (!(1-Epsilon<=sigma +alfa+beta && sigma+alfa+beta<=1+Epsilon)){
            continue;
        }
        if(!(tempT > 0)){
            continue;
        }        

        else{
            t=tempT;
            hit_flag = true;
            min_index = i;
        }

    }
    if(min_index != -1){
        result.normal_vector = triangles[min_index].indices.normal; 
    }
    
    if(hit_flag == true) {
        result.hit_time = t;
        Triangle triangle =  triangles[min_index];
        result.indices = triangle.indices;
        result.material_id = triangle.material_id;
        result.index = min_index;
    }
    return result;
}

float face_hit_time(const Face& face,const Vec3f& a, const Vec3f& b, const Vec3f& c, const Ray& ray) {
    const Vec3f &o = ray.o;
    const Vec3f &d = ray.d;
    float tempT=dotProduct(subtract(a,o),face.normal)/dotProduct(d,face.normal);
    if (tempT<0){
        return MAXFLOAT;
    }
    Vec3f hitPoint;
    hitPoint.x=o.x+d.x*tempT;
    hitPoint.y=o.y+d.y*tempT;
    hitPoint.z=o.z+d.z*tempT;
    Vec3f AC=subtract(c,a);
    float inv_AreaABC=1.0 / (magnitudeCalc(crossProduct(AC,face.edge_2)));
    Vec3f PB=subtract(b,hitPoint),PC=subtract(c,hitPoint);
    float alfa=magnitudeCalc(crossProduct(PB,PC))*inv_AreaABC;
    Vec3f PA=subtract(a,hitPoint);
    float beta=magnitudeCalc(crossProduct(PC,PA))*inv_AreaABC;
    float sigma=magnitudeCalc(crossProduct(PA,PB))*inv_AreaABC;
    
    if(!(-Epsilon<=alfa && alfa<=1+Epsilon)){
        return MAXFLOAT;    }
    if(!(-Epsilon<=beta && beta<=1+Epsilon)){
        return MAXFLOAT;
    }
    if(!(-Epsilon<=sigma && sigma<=1+Epsilon)){
        return MAXFLOAT;
    }
    if (!(1-Epsilon<=sigma +alfa+beta && sigma+alfa+beta<=1+Epsilon)){
        return MAXFLOAT;
    }
    if(!(tempT > 0)){
        return MAXFLOAT;
    }   
    else{
        return tempT;
    }    
}


Intersection mesh_hit(const Scene& scene,const Ray& ray){
    bool hit_flag = false;
    Intersection result;
    int min_mesh_index, min_face_index;
    result.type = MESH;
    result.hit_time = MAXFLOAT;
    Vec3f c,o=ray.o,d=ray.d;
    const vector<Mesh>& meshes = scene.meshes;
    const vector<Vec3f>& vertices = scene.vertex_data;
    for(int mesh_index = 0; mesh_index < meshes.size(); mesh_index++) {
        Mesh mesh = meshes[mesh_index];
        const vector<Face>& mesh_faces = mesh.faces;
        for (int face_index = 0; face_index < mesh_faces.size(); face_index++) {
            Face face = mesh_faces[face_index];
            int v0 = face.v0_id, v1 = face.v1_id, v2 = face.v2_id;
            Vec3f a = vertices[v0 - 1], b = vertices[v1 - 1], c = vertices[v2 - 1];
            float f_hit_time = face_hit_time(face,a, b, c, ray);
            if ( result.hit_time > f_hit_time ) {
                result.hit_time = f_hit_time;
                min_face_index = face_index;
                min_mesh_index = mesh_index;
                hit_flag = true;
            }
        }
    }
    if(hit_flag){
        Mesh mesh = meshes[min_mesh_index];
        Face face = mesh.faces[min_face_index];
        Vec3f a = vertices[face.v0_id - 1], b = vertices[face.v1_id - 1], c = vertices[face.v2_id - 1];
        Vec3f normal=normalizeVector(crossProduct(subtract(c,b),subtract(a,b)));
        result.normal_vector = normal;
        result.faces = mesh.faces;
        result.material_id = mesh.material_id;
        result.index = min_mesh_index;
    }
    return result;
}

Intersection rayIntersection(const Scene& scene,const Ray& ray){
    return firstIntersection(mesh_hit(scene, ray), firstIntersection(triangle_hit(scene, ray), sphere_hit(scene, ray)));
}


Vec3f calculateColorPixel(const Scene& scene, Ray ray,const Intersection& intersection, const Camera& currentCamera, int recursion_depth) {
    Vec3f pixelColor, intersection_point;
    float r, g, b;
    int material_id = intersection.material_id;
    int light_number = scene.point_lights.size();
    const vector<Material>& materials = scene.materials;
    intersection_point = add(ray.o, scalarMult(ray.d, intersection.hit_time)); 
    if(intersection.hit_time != MAXFLOAT) {
        r = scene.materials[material_id - 1].ambient.x * scene.ambient_light.x;
	    g = scene.materials[material_id - 1].ambient.y * scene.ambient_light.y;
	    b = scene.materials[material_id - 1].ambient.z * scene.ambient_light.z;
        for(int light_index = 0; light_index < light_number; light_index++){
            PointLight light = scene.point_lights[light_index]; // position & intensity
            Vec3f w_i = subtract(light.position, intersection_point);
            float distance_light_point = magnitudeCalc(w_i);
            w_i = normalizeVector(w_i);
            Vec3f wiEpsilon;
			wiEpsilon.x = w_i.x * scene.shadow_ray_epsilon;
			wiEpsilon.y = w_i.y * scene.shadow_ray_epsilon;
			wiEpsilon.z = w_i.z * scene.shadow_ray_epsilon;

            Ray shadow_ray;
            shadow_ray.o = add(intersection_point, wiEpsilon);
            shadow_ray.d = w_i;
            float hit_first_obj_time = rayIntersection(scene, shadow_ray).hit_time; 
            float light_hit = subtract(light.position, shadow_ray.o).x/shadow_ray.d.x;
            if(!(hit_first_obj_time < light_hit)){
                float cos_theta = max(float(0), dotProduct(w_i, intersection.normal_vector));
                Material material =  materials[intersection.material_id - 1];
                float r_sqr = distance_light_point*distance_light_point;
                r += (material.diffuse.x * cos_theta * light.intensity.x) / (r_sqr);// diffuse shading 
                g += (material.diffuse.y * cos_theta * light.intensity.y) / (r_sqr);// diffuse shading
                b += (material.diffuse.z * cos_theta * light.intensity.z) / (r_sqr);// diffuse shading
                Vec3f h = normalizeVector(subtract(w_i, ray.d));
                cos_theta = max(float(0), dotProduct(intersection.normal_vector, h));
                r += (material.specular.x * pow(cos_theta, material.phong_exponent) * light.intensity.x) / (r_sqr);// specular shading 
                g += (material.specular.y * pow(cos_theta, material.phong_exponent) * light.intensity.y) / (r_sqr);// specular shading
                b += (material.specular.z * pow(cos_theta, material.phong_exponent) * light.intensity.z) / (r_sqr);// specular shading
            }
            
        }

        // Mirrorness
        Vec3f refl_vector = {0,0,0};
        if(materials[material_id - 1].is_mirror && scene.max_recursion_depth > recursion_depth){
            float wi = -2 * dotProduct(ray.d, intersection.normal_vector);
            Vec3f wiNormal;
            wiNormal = add(scalarMult(intersection.normal_vector, wi), ray.d);
      	    wiNormal = normalizeVector(wiNormal);
      	    Vec3f wiEpsilon;
            wiEpsilon = scalarMult(wiNormal, scene.shadow_ray_epsilon);
            Ray refl;
            refl.o = add(intersection_point, wiEpsilon);
            refl.d = wiNormal;
            Intersection intersection_refl = rayIntersection(scene, refl);
            if(intersection_refl.hit_time != MAXFLOAT)
            {
                refl_vector = calculateColorPixel(scene, refl, intersection_refl, currentCamera, recursion_depth+1);
            }
      	    r += refl_vector.x * materials[material_id - 1].mirror.x;
      	    g += refl_vector.y * materials[material_id - 1].mirror.y;
      	    b += refl_vector.z * materials[material_id - 1].mirror.z;
   		
        }
    }else{
        if( !recursion_depth ){
            r = scene.background_color.x;
      	    g = scene.background_color.y;
      	    b = scene.background_color.z;
        }else{
            r = g = b = 0;
        }
    }
  	pixelColor.x = r;
  	pixelColor.y = g;
  	pixelColor.z = b;

    return pixelColor;
}






void* threadHandler(void* threadStruct){
    ThreadStruct *ts= (ThreadStruct*) threadStruct;
    int pixelNumber=0;
    
    for (int i=(*ts).beginI; i<(*ts).endI ; i++){
        
        for (int j=(*ts).beginJ; j<(*ts).endJ; j++){
            const Ray ray = generateRay(i, j, (*ts).currentCam);
            const Intersection ray_intersection = rayIntersection(scene, ray);
            pixelNumber=((*ts).w*i+j)*3;
            
            if(ray_intersection.hit_time == MAXFLOAT){
                image[pixelNumber] = scene.background_color.x;
                image[pixelNumber+1] = scene.background_color.y;
                image[pixelNumber+2] = scene.background_color.z;

            }else{
                Vec3f color_pixel = calculateColorPixel(scene, ray, ray_intersection, (*ts).currentCam, 0);
                color_pixel.x > 255 ? image[pixelNumber] = 255 : image[pixelNumber] = round(color_pixel.x);
                color_pixel.y > 255 ? image[pixelNumber + 1] = 255 : image[pixelNumber + 1] = round(color_pixel.y);
                color_pixel.z > 255 ? image[pixelNumber + 2] = 255 : image[pixelNumber + 2] = round(color_pixel.z);
            }
        }
    }
    return NULL;
}





int main(int argc, char* argv[])
{
    
    
    scene.loadFromXml(argv[1]);

    int cameraNum = scene.cameras.size();

    Camera currentCam;


    
    for (int camNum = 0; camNum < cameraNum; camNum++)
    {
        int pixelNumber = 0;
        currentCam = scene.cameras[camNum];
        w = currentCam.image_width;
        h = currentCam.image_height;


    
        image = new unsigned char[w * h * 3];
        int threadNum=4;
        int threadSize=h/4;
        ThreadStruct threadst[4];
        pthread_t threads[4];
        int t=0;
        int i=0;
        for (; i+threadSize<h; i+=threadSize,t++){
            ThreadStruct ts;
            ts.currentCam=currentCam;
            
            ts.h=h;
            ts.w=w;
            ts.beginI=i;
            ts.beginJ=0;
            ts.endI=i+threadSize;
            ts.endJ=w;
            threadst[t]=ts;
        }
        ThreadStruct ts;
        ts.currentCam=currentCam;
        ts.h=h;
        ts.w=w;
        ts.beginI=i;
        ts.beginJ=0;
        ts.endI=h;
        ts.endJ=w;
        threadst[t]=ts;
        for (int i=0; i<threadNum;i++){
            pthread_create(&threads[i],NULL,&threadHandler,(void*)&threadst[i]);
        }

        for (int i=0; i<threadNum;i++){
            
            pthread_join(threads[i],NULL);
        }
        write_ppm(currentCam.image_name.c_str(), image, w, h);
        free(image);
    }
    
}
