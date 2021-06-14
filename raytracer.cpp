#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "model.h"
#include "geometry.h"

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;

int surface_width, surface_height;
std::vector<Vec3f> surface;
Model ico("../icosahedron.obj");

Vec3f ico_center(0, 0, -16);

struct Light {
    Light(const Vec3f &p, const float i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Material {
    Material(const float r, const Vec4f &a, const Vec3f &color, const float spec, bool cr = false, bool sr = false) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec), crystal(cr), surface(sr) {}
    Material() : refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent(), crystal(false) {}
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
    bool crystal;
    bool surface;
};

Material mat;

struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float r, const Material &m) : center(c), radius(r), material(m) {};

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

struct Cube {
    Vec3f center;
    Vec3f sizes;
    float x, y, z;
    Material material;
    Cube(const Vec3f &c, const Vec3f &s, const Material &m): center(c), sizes(s), material(m) {};

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0, Vec3f& normal, float &u, float &v) const {
        float x = sizes.x;
        float y = sizes.y;
        float z = sizes.z;
        float l = std::numeric_limits<float>::max();
        Vec3f worldPos;
        float t;
        t = (-z / 2.0 + center.z - orig.z) / dir.z;
        worldPos = orig + dir * t - center;
        if (l > t && t > 0.0 && std::abs(worldPos.x) < x/2. && std::abs(worldPos.y) < y/2.) {
            l = t;
            t0 = l;
            u = (worldPos.x + x/2) / x;
            v = (worldPos.y + y/2) / y;
            normal = Vec3f(0, 0, -1);
        }
        t = (z / 2.0 + center.z - orig.z) / dir.z;
        worldPos = orig + dir * t - center;
        if (l > t && t > 0.0 && std::abs(worldPos.x) < x/2. && std::abs(worldPos.y) < y/2.) {
            l = t;
            t0 = l;
            u = (worldPos.x + x/2) / x;
            v = (worldPos.y + y/2) / y;
            normal = Vec3f(0, 0, 1);
        }
        t = (-x / 2.0 + center.x - orig.x) / dir.x;
        worldPos = orig + dir * t - center;
        if (l > t && t > 0.0 && std::abs(worldPos.y) < y/2. && std::abs(worldPos.z) < z/2.) {
            l = t;
            t0 = l;
            u = (worldPos.z + z/2) / z;
            v = (worldPos.y + y/2) / y;
            normal = Vec3f(-1, 0, 0);
        }
        t = (x / 2.0 + center.x - orig.x) / dir.x;
        worldPos = orig + dir * t - center;
        if (l > t && t > 0.0 && std::abs(worldPos.y) < y/2. && std::abs(worldPos.z) < z/2.) {
            l = t;
            t0 = l;
            u = (worldPos.z + z/2) / z;
            v = (worldPos.y + y/2) / y;
            normal = Vec3f(1, 0, 0);
        }
        
        t = (-y / 2.0 + center.y - orig.y) / dir.y;
        worldPos = orig + dir * t - center;
        if (l > t && t > 0.0 && std::abs(worldPos.x) < x/2. && std::abs(worldPos.z) < z/2.) {
            l = t;
            t0 = l;
            u = (worldPos.z + z/2) / z;
            v = (worldPos.x + x/2) / x;
            normal = Vec3f(0, -1, 0);
        }
        t = (y / 2.0 + center.y - orig.y) / dir.y;
        worldPos = orig + dir * t - center;
        if (l > t && t > 0.0 && std::abs(worldPos.x) < x/2. && std::abs(worldPos.z) < z/2.) {
            l = t;
            t0 = l;
            u = (worldPos.z + z/2) / z;
            v = (worldPos.x + x/2) / x;
            normal = Vec3f(0, 1, 0);
        }
        return l < 1000;
    }
};

std::vector<Sphere> spheres;
std::vector<Light>  lights;
std::vector<Cube> cubes;


Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t, const float eta_i=1.f) { 
    float cosi = - std::max(-1.f, std::min(1.f, I*N));
    if (cosi<0) return refract(I, -N, eta_i, eta_t);
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k<0 ? Vec3f(1,0,0) : I*eta + N*(eta*cosi - sqrtf(k));
}

void fresnel(const Vec3f &I, const Vec3f &N, const float &ior, float &kr) 
{ 
    float cosi = std::max(-1.f, std::min(1.f, I*N));
    float etai = 1, etat = ior; 
    if (cosi > 0) { std::swap(etai, etat); } 
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    if (sint >= 1) { 
        kr = 1; 
    } else { 
        float cost = sqrtf(std::max(0.f, 1 - sint * sint)); 
        cosi = fabsf(cosi); 
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
        kr = (Rs * Rs + Rp * Rp) / 2; 
    }
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N, Material &material) {

    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }
    for (size_t i=0; i < (size_t)ico.nfaces(); i++) {
        float dist_i;
        Vec3f norm;
        float u, v;
        if (ico.ray_triangle_intersect(i, orig, dir, dist_i, norm, u, v) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = norm.normalize();
            material = mat;
            material.crystal = true;
        }
    }
    for (size_t i=0; i < cubes.size(); i++) {
        float dist_i;
        Vec3f norm;
        float u, v;
        if (cubes[i].ray_intersect(orig, dir, dist_i, norm, u, v) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = norm.normalize();
            material = cubes[i].material;
            if (material.surface &&
            (int)(u * surface_width) + (int)(v * surface_height) * surface_width >= 0 &&
            (int)(u * surface_width) + (int)(v * surface_height) * surface_width < surface_height * surface_width) {
                material.diffuse_color = surface[(int)(u * surface_width) + (int)(v * surface_height) * surface_width];
            }
        }
    }
    return spheres_dist < 1000;
}

bool light(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N, Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }
    for (size_t i=0; i < cubes.size(); i++) {
        float dist_i;
        Vec3f norm;
        float u, v;
        if (cubes[i].ray_intersect(orig, dir, dist_i, norm, u, v) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = norm.normalize();
            material = cubes[i].material;
        }
    }
    return spheres_dist<1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, size_t depth=0) {
    Vec3f point, N;
    Material material;

    if (depth>3 || !scene_intersect(orig, dir, point, N, material)) {
        Sphere env(Vec3f(0,0,0), 100, Material());
        float dist = 0;
        env.ray_intersect(orig, dir, dist);
        Vec3f p = orig+dir*dist;
        int a = (atan2(p.z, p.x)/(2*M_PI) + .5)*envmap_width;
        int b = acos(p.y/100)/M_PI*envmap_height;
        return envmap[a+b*envmap_width];
    }

    if (material.crystal) {
        Vec3f refractionColor = Vec3f(0, 0, 0);
        float q;
        fresnel(dir, N, material.refractive_index, q);
        bool outside = dir * N < 0;
        Vec3f bias = N * 1e-3;
        if (q < 1) {
            Vec3f refractionDirection = refract(dir, N, material.refractive_index).normalize(); 
            Vec3f refractionRayOrig = outside ? point - bias : point + bias; 
            refractionColor = cast_ray(refractionRayOrig, refractionDirection, depth + 1); 
        }
        Vec3f reflectionDirection = reflect(dir, N).normalize(); 
        Vec3f reflectionRayOrig = outside ? point + bias : point - bias; 
        Vec3f reflectionColor = cast_ray(reflectionRayOrig, reflectionDirection, depth + 1); 

        float diffuse_light_intensity = 0, specular_light_intensity = 0;
        for (size_t i=0; i<lights.size(); i++) {
            Vec3f light_dir      = (lights[i].position - point).normalize();
            float light_distance = (lights[i].position - point).norm();

            Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
            Vec3f shadow_pt, shadow_N;
            Material tmpmaterial;
            if (light(shadow_orig, light_dir, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
                continue;

            specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
        }
        return Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflectionColor * q + refractionColor * (1 - q);
    }
    
    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, depth + 1);
    Vec3f refract_color = cast_ray(refract_orig, refract_dir, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (light(shadow_orig, light_dir, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflect_color *material.albedo[2] + refract_color*material.albedo[3];
}

Vec3f rotation(const Vec3f& v, const Vec3f& k, double theta)
{
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    Vec3f rotated = (v * cos_theta) + (cross(k, v) * sin_theta) + (k * (k * v)) * (1 - cos_theta);

    return rotated;
}

void render(int width, int height) {
    const float fov      = M_PI/3.;
    std::vector<Vec3f> framebuffer(width*height);

    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float dir_x =  (i + 0.5) -  width/2.;
            float dir_y = -(j + 0.5) + height/2. ;
            float dir_z = -height/(2.*tan(fov/2.));
            Vec3f dir(dir_x, dir_y, dir_z);
            
            dir = rotation(dir, Vec3f(0, 1, 0).normalize(), -M_PI / 20);
            dir = rotation(dir, Vec3f(1, 0, 0).normalize(), -M_PI / 20);

            Vec3f dir0(dir.x + 0.5, dir.y, dir.z);
            Vec3f dir1(dir.x, dir.y + 0.5, dir.z);
            Vec3f dir2(dir.x + 0.5, dir.y + 0.5, dir.z);
            Vec3f pos(-3, 2, 0);

            framebuffer[i+j*width] = (cast_ray(pos, dir.normalize()) + 
            cast_ray(pos, dir0.normalize()) + cast_ray(pos, dir1.normalize()) +
            cast_ray(pos, dir2.normalize())) * 0.25;
        }
    }

    std::vector<unsigned char> pixmap(width*height*3);
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            pixmap[i*3+j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    stbi_write_jpg("316_Klimenko_v5v3.jpg", width, height, 3, pixmap.data(), 100);
}

int main(int argc, char *argv[]) {
    int n = -1;
    unsigned char *pixmap = stbi_load("../envmap.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    envmap = std::vector<Vec3f>(envmap_width*envmap_height);
    for (int j = envmap_height-1; j>=0 ; j--) {
        for (int i = 0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);
    pixmap = stbi_load("../surface.jpg", &surface_width, &surface_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the surface texture" << std::endl;
        return -1;
    }
    surface = std::vector<Vec3f>(surface_width*surface_height);
    for (int j = surface_height-1; j>=0 ; j--) {
        for (int i = 0; i<surface_width; i++) {
            surface[i+j*surface_width] = Vec3f(pixmap[(i+j*surface_width)*3+0], pixmap[(i+j*surface_width)*3+1], pixmap[(i+j*surface_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);

    Material glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125.);
    Material littlemsu(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.3, 0.2),   10.);
    Material mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);
    Material msu_mat(1.0, Vec4f(0.5,  0.3, 0.0, 0.5), Vec3f(0.4, 0.4, 0.3),   50.);
    Material surface_mat(1.0, Vec4f(0.4,  0.3, 0.0, 0.0), Vec3f(0.4, 0.4, 0.3),   50., false, true);

    mat = glass;


    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));

    //////////////////////////////////////////////////////////////

    float mult = 4;
    Vec3f offset(0, 0, -16);

    for (int i = 0; i < ico.nverts(); ++i) {
        ico.point(i) = ico.point(i) * mult;
        ico.point(i) = ico.point(i) + offset;
    }

    /////////////////////////////////////////////

    Material building = msu_mat;
    Vec3f place = Vec3f (0, 0, -16);
    float build_mult = 0.3;

    Vec3f a(0, 0, 0), as(3, 4, 3); // башня
    Vec3f f(0, 3, 0), fs(2, 2, 2);
    Vec3f g(0, 4.625, 0), gs(1, 1.25, 1);

    Vec3f b(-2.75, -0.75, 0), bs(2.5, 2.5, 1); // крылья
    Vec3f c(2.75, -0.75, 0), cs(2.5, 2.5, 1);

    Vec3f d(-4.75, -0.5, 0), ds(1.5, 3, 4); // боковые
    Vec3f e(4.75, -0.5, 0), es(1.5, 3, 4);

    Vec3f h(0, 5.875, 0), hs(0.125, 1.25, 0.125); // шпиль
    spheres.push_back(Sphere(Vec3f(0, 6.65, 0) * build_mult + place, 0.2 * build_mult, building));

    Vec3f i(-4.75, 1.5, 1.25), is(1, 1, 1);
    Vec3f j(-4.75, 1.5, -1.25), js(1, 1, 1);
    Vec3f k(4.75, 1.5, 1.25), ks(1, 1, 1);
    Vec3f l(4.75, 1.5, -1.25), ls(1, 1, 1);

    Vec3f m(-4.75, 2.25, 1.25), ms(0.5, 0.5, 0.5);
    Vec3f q(-4.75, 2.25, -1.25), qs(0.5, 0.5, 0.5);
    Vec3f o(4.75, 2.25, 1.25), os(0.5, 0.5, 0.5);
    Vec3f p(4.75, 2.25, -1.25), ps(0.5, 0.5, 0.5);

    cubes.push_back(Cube(a * build_mult + place, as * build_mult, building));
    cubes.push_back(Cube(b * build_mult + place, bs * build_mult, building));
    cubes.push_back(Cube(c * build_mult + place, cs * build_mult, building));
    cubes.push_back(Cube(d * build_mult + place, ds * build_mult, building));
    cubes.push_back(Cube(e * build_mult + place, es * build_mult, building));
    cubes.push_back(Cube(f * build_mult + place, fs * build_mult, building));
    cubes.push_back(Cube(g * build_mult + place, gs * build_mult, building));
    cubes.push_back(Cube(h * build_mult + place, hs * build_mult, building));
    cubes.push_back(Cube(i * build_mult + place, is * build_mult, building));
    cubes.push_back(Cube(j * build_mult + place, js * build_mult, building));
    cubes.push_back(Cube(k * build_mult + place, ks * build_mult, building));
    cubes.push_back(Cube(l * build_mult + place, ls * build_mult, building));
    cubes.push_back(Cube(m * build_mult + place, ms * build_mult, building));
    cubes.push_back(Cube(q * build_mult + place, qs * build_mult, building));
    cubes.push_back(Cube(o * build_mult + place, os * build_mult, building));
    cubes.push_back(Cube(p * build_mult + place, ps * build_mult, building));

    float build_mult2 = 0.2;
    Vec3f place2(3, -3.4, -11);
    Material building2 = littlemsu;

    cubes.push_back(Cube(a * build_mult2 + place2, as * build_mult2, building2));
    cubes.push_back(Cube(b * build_mult2 + place2, bs * build_mult2, building2));
    cubes.push_back(Cube(c * build_mult2 + place2, cs * build_mult2, building2));
    cubes.push_back(Cube(d * build_mult2 + place2, ds * build_mult2, building2));
    cubes.push_back(Cube(e * build_mult2 + place2, es * build_mult2, building2));
    cubes.push_back(Cube(f * build_mult2 + place2, fs * build_mult2, building2));
    cubes.push_back(Cube(g * build_mult2 + place2, gs * build_mult2, building2));
    cubes.push_back(Cube(h * build_mult2 + place2, hs * build_mult2, building2));
    cubes.push_back(Cube(i * build_mult2 + place2, is * build_mult2, building2));
    cubes.push_back(Cube(j * build_mult2 + place2, js * build_mult2, building2));
    cubes.push_back(Cube(k * build_mult2 + place2, ks * build_mult2, building2));
    cubes.push_back(Cube(l * build_mult2 + place2, ls * build_mult2, building2));
    cubes.push_back(Cube(m * build_mult2 + place2, ms * build_mult2, building2));
    cubes.push_back(Cube(q * build_mult2 + place2, qs * build_mult2, building2));
    cubes.push_back(Cube(o * build_mult2 + place2, os * build_mult2, building2));
    cubes.push_back(Cube(p * build_mult2 + place2, ps * build_mult2, building2));
    spheres.push_back(Sphere(Vec3f(0, 6.65, 0) * build_mult2 + place2, 0.2 * build_mult2, building)); // "звезда"


    // surface
    cubes.push_back(Cube(Vec3f(0, -4, -16), Vec3f(20, 0.5, 20), surface_mat));

    cubes.push_back(Cube(Vec3f(0, -3.5, -16), Vec3f(4.5, 0.5, 4.5), surface_mat));

    spheres.push_back(Sphere(Vec3f(-5, -2, -11), 1.75, mirror));

    int width = 512, height = 512;
    for (int i = 1; i < argc - 1; ++i) {
        if (argv[i][0] == '-' && argv[i][1] == 'w') {
            width = height = std::stoi(argv[i + 1]);
        }
    }

    render(width, height);

    return 0;
}

