#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include "model.h"

// fills verts and faces arrays, supposes .obj file to have "f " entries without slashes
Model::Model(const char *filename) : verts(), faces() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v[i];
            verts.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            Vec3i f;
            int idx, cnt=0;
            iss >> trash;
            while (iss >> idx) {
                idx--; // in wavefront obj all indices start at 1, not zero
                f[cnt++] = idx;
            }
            if (3==cnt) faces.push_back(f);
        }
    }
    std::cerr << "# v# " << verts.size() << " f# "  << faces.size() << std::endl;

    Vec3f min, max;
    get_bbox(min, max);
}

// Moller and Trumbore
bool Model::ray_triangle_intersect(const int &fi, const Vec3f &orig, const Vec3f &dir, float &tnear, Vec3f &N, float &u, float &v) {
    Vec3f v0v1 = point(vert(fi,1)) - point(vert(fi,0));
    Vec3f v0v2 = point(vert(fi,2)) - point(vert(fi,0));
    Vec3f pvec = cross(dir, v0v2); 
    float det = v0v1 * pvec; 
    N = cross(v0v1, v0v2);
    
    if (fabs(det) < 1e-6) return false; 
    float invDet = 1 / det; 

    Vec3f tvec = orig - point(vert(fi,0));
    u = (tvec * pvec) * invDet; 
    if (u < 0 || u > 1) return false; 

    Vec3f qvec = cross(tvec, v0v1); 
    v = (dir * qvec) * invDet; 
    if (v < 0 || u + v > 1) return false; 

    tnear = (v0v2 * qvec) * invDet;

    return tnear > 0; 
}


int Model::nverts() const {
    return (int)verts.size();
}

int Model::nfaces() const {
    return (int)faces.size();
}

void Model::get_bbox(Vec3f &min, Vec3f &max) {
    min = max = verts[0];
    for (int i=1; i<(int)verts.size(); ++i) {
        for (int j=0; j<3; j++) {
            min[j] = std::min(min[j], verts[i][j]);
            max[j] = std::max(max[j], verts[i][j]);
        }
    }
    std::cerr << "bbox: [" << min << " : " << max << "]" << std::endl;
}

const Vec3f &Model::point(int i) const {
    assert(i>=0 && i<nverts());
    return verts[i];
}

Vec3f &Model::point(int i) {
    assert(i>=0 && i<nverts());
    return verts[i];
}

int Model::vert(int fi, int li) const {
    assert(fi>=0 && fi<nfaces() && li>=0 && li<3);
    return faces[fi][li];
}

std::ostream& operator<<(std::ostream& out, Model &m) {
    for (int i=0; i<m.nverts(); i++) {
        out << "v " << m.point(i) << std::endl;
    }
    for (int i=0; i<m.nfaces(); i++) {
        out << "f ";
        for (int k=0; k<3; k++) {
            out << (m.vert(i,k)+1) << " ";
        }
        out << std::endl;
    }
    return out;
}

