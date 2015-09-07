#include "scene.h"

// intersection record
struct intersection3f {
    bool        hit;        // whether it hits something
    float       ray_t;      // ray parameter for the hit
    vec3f       pos;        // hit position
    vec3f       norm;       // hit normal
    Material*   mat;        // hit material

    // constructor (defaults to no intersection)
    intersection3f() : hit(false) { }
};

#define ray3f_epsilon 0.0005f
#define ray3f_rayinf 1000000.0f

// 3D Ray
struct ray3f {
    vec3f e;        // origin
    vec3f d;        // direction
    float tmin;     // min t value
    float tmax;     // max t value

    // Default constructor
    ray3f() : e(zero3f), d(z3f), tmin(ray3f_epsilon), tmax(ray3f_rayinf) { }

    // Element-wise constructor
    ray3f(const vec3f& e, const vec3f& d) :
        e(e), d(d), tmin(ray3f_epsilon), tmax(ray3f_rayinf) { }

    // Element-wise constructor
    ray3f(const vec3f& e, const vec3f& d, float tmin, float tmax) :
        e(e), d(d), tmin(tmin), tmax(tmax) { }

    // Eval ray at a specific t
    vec3f eval(float t) const { return e + d * t; }

    // Create a ray from a segment
    static ray3f make_segment(const vec3f& a, const vec3f& b) { return ray3f(a,normalize(b-a),ray3f_epsilon,dist(a,b)-2*ray3f_epsilon); }
};

// transform a ray by a frame
inline ray3f transform_ray(const frame3f& f, const ray3f& v) {
    return ray3f(transform_point(f,v.e), transform_vector(f,v.d), v.tmin, v.tmax);
}
// transform a ray by a frame inverse
inline ray3f transform_ray_inverse(const frame3f& f, const ray3f& v) {
    return ray3f(transform_point_inverse(f,v.e),transform_vector_inverse(f,v.d),v.tmin,v.tmax);
}


// intersects the scene and return the first intrerseciton
intersection3f intersect(Scene* scene, ray3f ray) {

    // create a default intersection record to be returned
    auto intersection = intersection3f();

    //set the closest t value to be just farther than the t max
    float closest_t = ray.tmax +1;
    intersection3f closest_int = intersection3f();

    //foreach surface
    for(int i=0; i < scene->surfaces.size(); i++){
        auto cur_shape = scene->surfaces[i];
        // if it is a quad

        if(cur_shape->isquad == true){
            //compute ray intersection (and ray parameter), continue if not hit
            vec3f n = cur_shape->frame.z;
            vec3f C = cur_shape->frame.o;

            vec3f E = ray.e;
            float R = cur_shape->radius;

            if(dot(normalize(ray.d), n) == 0){continue;}

            float t = dot((C-E), n)/(dot(normalize(ray.d), n));

            if(t < ray.tmin || t >ray.tmax){continue;}

            vec3f PT = ray.eval(t);         //find intersection point with quad
            vec3f localPT = transform_point_inverse(cur_shape->frame, PT);      //transform intersection point to shape frame

            if((abs(localPT.x) - cur_shape->frame.o.x) < R && (abs(localPT.y) -cur_shape->frame.o.y) < R){        //check to see if local point is within dimentions
                if(t < closest_t){
                    closest_t = t;
                    closest_int.hit = true;
                    closest_int.mat = cur_shape->mat;
                    closest_int.ray_t = t;
                    closest_int.pos = PT;           //return PT
                    closest_int.norm = normalize(n);
                }
            }


        }

        // if it is a sphere
        else{
            // compute ray intersection (and ray parameter), continue if not hit
            vec3f C = cur_shape->frame.o; //sphere center
            vec3f E = ray.e; //camera origin
            float R = cur_shape->radius; //sphere radius


            float a = lengthSqr(normalize(ray.d));
            float b = dot(2*normalize(ray.d), (E-C));
            float c = lengthSqr(E-C)-sqr(R);

            float det = sqr(b) - 4*a*c;


            if(det < 0){continue;}
            //calculate t1 and t2
            float t1 = (-b + sqrt(det))/(2*a);
            float t2 = (-b - sqrt(det))/(2*a);
            //             just grab only the first hit
            float t = fminf(t1, t2);
            //check if computed param is within ray.tmin and ray.tmax
            if(t >= ray.tmin && t <= ray.tmax){
                if(t < closest_t){
                    closest_t = t;
                    closest_int.hit = true;
                    closest_int.mat = cur_shape->mat;
                    closest_int.ray_t = t;
                    closest_int.pos = ray.eval(t);
                    closest_int.norm = normalize(closest_int.pos - C);
                }
            }



        }
    }
    intersection = closest_int;
    return intersection;
}

// compute the color corresponing to a ray by raytracing
vec3f raytrace_ray(Scene* scene, ray3f ray) {

    // get scene intersection
    intersection3f myint = intersect(scene, ray);
    // if not hit, return background
    if(myint.hit == false){

        return scene->background;

    }

    vec3f kd = myint.mat->kd; //material response

    // accumulate color starting with ambient
    vec3f C = scene->ambient * kd; //*kd


    // foreach light
    for(auto &light : scene->lights){
        // compute light response
        vec3f ks = myint.mat->ks;
        // compute light direction
        vec3f ld = normalize(light->frame.o - myint.pos);
        // compute the material response (brdf*cos) ABOVE

        // check for shadows and accumulate if needed
        ray3f shaddow_ray = ray3f(light->frame.o, -ld);     //-light direction is the direction of shaddow
        shaddow_ray.tmax = length(light->frame.o - myint.pos) - ray3f_epsilon;      //-ray3f_epsilon to
        intersection3f shaddow_int = intersect(scene, shaddow_ray);

        if(shaddow_int.hit){continue;}

        vec3f h = normalize(ld - normalize(ray.d));
        vec3f n = myint.norm;
        float nexp = myint.mat->n;

        // return the accumulated color (for now zero)
        //using blinn phong:
        vec3f Li= light->intensity/(lengthSqr(light->frame.o - myint.pos));
        vec3f spec = (kd + ks*pow(fmax(0, dot(n, h)), nexp))*abs(dot(n, ld));


        C += Li*spec;
    }
    // if the material has reflections
    vec3f AR;
    if(myint.mat->kr != zero3f){
        // create the reflection ray
        ray3f ref_ray = ray3f(myint.pos, normalize(ray.d)+ 2*(dot(-1*normalize(ray.d), myint.norm))*myint.norm);
        // accumulate the reflected light (recursive call) scaled by the material reflection

        AR = myint.mat->kr * raytrace_ray(scene, ref_ray);

    }

    // return the accumulated color (for now zero)


    return C + AR;
}

// raytrace an image
image3f raytrace(Scene* scene) {
    // allocate an image of the proper size
    auto image = image3f(scene->image_width, scene->image_height);


    int width_pixels = image.width();
    int height_pixels = image.height();

    // if no anti-aliasing
    if(scene->image_samples == 1){
        // foreach image row (go over image height)
        for(float p_h = 0; p_h <height_pixels; p_h++){
            // foreach pixel in the row (go over image width)
            for(float p_w = 0; p_w < width_pixels; p_w++){

                // compute ray-camera parameters (u,v) for the pixel
                float v = (p_h +0.5)/height_pixels;
                float u = (p_w + 0.5)/width_pixels;

                // compute camera ray
                Camera* cam = scene->camera;
                frame3f frame = cam->frame;
                vec3f Q = (u-0.5)*cam->width*frame.x+ (v-0.5)*cam->height*frame.y -cam->dist*frame.z;
                vec3f cam_origin = frame.o;
                ray3f curr_ray =ray3f(cam_origin, normalize(Q));
                // set pixel to the color raytraced with the ray
                image.at(p_w, p_h)=raytrace_ray(scene, curr_ray);

            }
        }
    }

    else{
        // foreach image row (go over image height)
        for(float p_h = 0; p_h <height_pixels; p_h++){
            // foreach pixel in the row (go over image width)
            for(float p_w = 0; p_w < width_pixels; p_w++){

                // init accumulated color
                vec3f C = zero3f;
                // foreach sample in y
                for(float p_hh = 0; p_hh < scene->image_samples; p_hh++){
                    // foreach sample in x
                    for(float p_ww = 0; p_ww < scene->image_samples; p_ww++){
                        // compute ray-camera parameters (u,v) for the pixel and the sample
                        // compute camera ray
                        // set pixel to the color raytraced with the ray

                        float v = (p_h +(p_hh + 0.5)/scene->image_samples)/height_pixels;
                        float u = (p_w + (p_ww + 0.5)/scene->image_samples)/width_pixels;

                        // compute camera ray
                        Camera* cam = scene->camera;
                        frame3f frame = cam->frame;
                        vec3f Q = (u-0.5)*cam->width*frame.x+ (v-0.5)*cam->height*frame.y -cam->dist*frame.z;
                        vec3f cam_origin = frame.o;
                        ray3f curr_ray =ray3f(cam_origin, normalize(Q));
                        // set pixel to the color raytraced with the ray
                        C += raytrace_ray(scene, curr_ray);
                    }
                }
                image.at(p_w, p_h) = C/(scene->image_samples * scene->image_samples);
            }
        }
    }
    // done
    return image;
}

// runs the raytrace over all tests and saves the corresponding images
int main(int argc, char** argv) {
    auto args = parse_cmdline(argc, argv,
    { "01_raytrace", "raytrace a scene",
      {  {"resolution",     "r", "image resolution", typeid(int),    true,  jsonvalue()}  },
      {  {"scene_filename", "",  "scene filename",   typeid(string), false, jsonvalue("scene.json")},
         {"image_filename", "",  "image filename",   typeid(string), true,  jsonvalue("")}  }
                              });

    // generate/load scene either by creating a test scene or loading from json file
    string scene_filename = args.object_element("scene_filename").as_string();
    Scene *scene = nullptr;
    if(scene_filename.length() > 9 and scene_filename.substr(0,9) == "testscene") {
        int scene_type = atoi(scene_filename.substr(9).c_str());
        scene = create_test_scene(scene_type);
        scene_filename = scene_filename + ".json";
    } else {
        scene = load_json_scene(scene_filename);
    }
    error_if_not(scene, "scene is nullptr");

    auto image_filename = (args.object_element("image_filename").as_string() != "") ?
                args.object_element("image_filename").as_string() :
                scene_filename.substr(0,scene_filename.size()-5)+".png";

    if(not args.object_element("resolution").is_null()) {
        scene->image_height = args.object_element("resolution").as_int();
        scene->image_width = scene->camera->width * scene->image_height / scene->camera->height;
    }

    message("rendering %s...\n", scene_filename.c_str());
    auto image = raytrace(scene);

    message("writing to png...\n");
    write_png(image_filename, image, true);

    delete scene;
    message("done\n");
}
