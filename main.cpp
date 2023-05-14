#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"
#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <list>

std::random_device rd;
thread_local std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }

    double norm() const {
        return sqrt(norm2());
    }

    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }

    double operator[](int i) const { return data[i]; };

    double &operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector &a, const Vector &b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector &a, const Vector &b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(const double a, const Vector &b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector &a, const double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

Vector operator*(const Vector &a, const Vector &b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector &a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}

Vector operator+=(Vector &a, const Vector &b) {
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
    return a;
}

double dot(const Vector &a, const Vector &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector &a, const Vector &b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
    Ray(const Vector &O, const Vector &u) : origin(O), direction(u) {};

    Vector origin;
    Vector direction;
};

class Intersection {
public:
    Intersection(Vector P = Vector(0, 0, 0), Vector N = Vector(0, 0, 0), double t = std::numeric_limits<double>::max(), Vector rho = Vector(0, 0, 0))
        : P(P), N(N), t(t), rho(rho) {};

    Vector P;
    Vector N;
    double t;
    Vector rho;
    bool isMirror;
};

class BoundingBox {
    public:

        BoundingBox() {
            m = Vector(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
            M = Vector(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
        }

        bool intersect(const Ray& r) const {
            double tx1 = (m[0] - r.origin[0]) / r.direction[0];
            double tx2 = (M[0] - r.origin[0]) / r.direction[0];
            double txmin = std::min(tx1, tx2);
            double txmax = std::max(tx1, tx2);

            double ty1 = (m[1] - r.origin[1]) / r.direction[1];
            double ty2 = (M[1] - r.origin[1]) / r.direction[1];
            double tymin = std::min(ty1, ty2);
            double tymax = std::max(ty1, ty2);

            double tz1 = (m[2] - r.origin[2]) / r.direction[2];
            double tz2 = (M[2] - r.origin[2]) / r.direction[2];
            double tzmin = std::min(tz1, tz2);
            double tzmax = std::max(tz1, tz2);

            double actualT = std::max(txmin, std::max(tymin, tzmin));
            double furthestT = std::min(txmax, std::min(tymax, tzmax));
            if (furthestT > 0 && furthestT > actualT) return true;
            return false;
        }

        void expand(const Vector& v) {
            for (int i = 0; i < 3; i++) {
                m[i] = std::min(m[i], v[i]);
                M[i] = std::max(M[i], v[i]);
            }
        }

        int longestAxis() const {
            Vector d = M - m; // Extent of the bounding box along each axis
            if (d[0] > d[1] && d[0] > d[2])
                return 0;
            else if (d[1] > d[2])
                return 1;
            else
                return 2;
        }

        Vector m, M;
};

class Object{
    public:
        Object(Vector rho, bool mirror) : rho(rho), mirror(mirror) {};
        virtual bool intersect(const Ray &r, Intersection &I) const = 0;
        Vector rho;
        bool mirror;
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class BVHNode {
public:
    BoundingBox bounding_box;
    BVHNode* left;
    BVHNode* right;
    size_t start;
    size_t end;

    BVHNode() : left(nullptr), right(nullptr), start(-1), end(-1) {}

    ~BVHNode() {
        delete left;
        delete right;
    }
};

class TriangleMesh: public Object {
public:

    TriangleMesh(Vector rho, bool mirror) : Object(rho, mirror), rho(rho), mirror(mirror) {};

    ~TriangleMesh() {
        delete root;
    }

    bool intersect_triangle(const Ray& ray, size_t i, Intersection &intersection) const {
        const TriangleIndices &tri = indices[i];
        const Vector &v0 = vertices[tri.vtxi];
        const Vector &v1 = vertices[tri.vtxj];
        const Vector &v2 = vertices[tri.vtxk];

        Vector edge1 = v1 - v0;
        Vector edge2 = v2 - v0;
        Vector h = cross(ray.direction,edge2);
        double a = dot(edge1,h);

        if (std::abs(a) < 1e-8) {
            return false;
        }

        double f = 1.0 / a;
        Vector s = ray.origin - v0;
        double u = f * dot(s,h);

        if (u < 0.0 || u > 1.0) {
            return false;
        }

        Vector q = cross(s, edge1);
        double v = f * dot(ray.direction,q);

        if (v < 0.0 || u + v > 1.0) {
            return false;
        }

        double t = f * dot(edge2,q);

        if (t > 1e-8 && t < intersection.t) {
            intersection.t = t;
            intersection.P = ray.origin + t * ray.direction;
            intersection.N = cross(edge1, edge2);
            intersection.N.normalize();
            return true;
        }

        return false;
    }

    BVHNode* build_bvh(size_t start, size_t end, int min_triangles) {
        BVHNode* node = new BVHNode();
        compute_bounding_box(node->bounding_box, start, end);

        size_t num_triangles = end - start;
        if (num_triangles <= min_triangles) {
            node->start = start;
            node->end = end;
        } else {
            int axis = node->bounding_box.longestAxis();
            auto compare_triangles = [&](const TriangleIndices &a, const TriangleIndices &b) {
                return vertices[a.vtxi][axis] + vertices[a.vtxj][axis] + vertices[a.vtxk][axis] <
                    vertices[b.vtxi][axis] + vertices[b.vtxj][axis] + vertices[b.vtxk][axis];
            };
            size_t mid = start + (end - start) / 2;
            std::nth_element(indices.begin() + start, indices.begin() + mid, indices.begin() + end, compare_triangles);

            node->left = build_bvh(start, mid, min_triangles);
            node->right = build_bvh(mid, end, min_triangles);
        }

        return node;
    }

    void intersectBVH(const Ray &ray, BVHNode *node, Intersection &intersection) const {
        if (!node || !node->bounding_box.intersect(ray)) {
            return;
        }
        if (node->left == nullptr && node->right == nullptr) {
            for (size_t i = node->start; i < node->end; ++i) {
                Intersection temp_intersection = intersection;
                if (intersect_triangle(ray, i, temp_intersection) && temp_intersection.t < intersection.t) {
                    intersection = temp_intersection;
                }
            }
        } else {
            intersectBVH(ray, node->left, intersection);
            intersectBVH(ray, node->right, intersection);
        }
    }

    bool intersect(const Ray& ray, Intersection &intersection) const {
        intersectBVH(ray, root, intersection);

        if (intersection.t < std::numeric_limits<double>::max()) {
            intersection.rho = rho;
            intersection.isMirror = mirror;
            return true;
        } else {
            return false;
        }
    }

    void compute_bounding_box(BoundingBox &bbox, size_t start, size_t end) {
        for (size_t i = start; i < end; i++) {
            const TriangleIndices &tri = indices[i];
            bbox.expand(vertices[tri.vtxi]);
            bbox.expand(vertices[tri.vtxj]);
            bbox.expand(vertices[tri.vtxk]);
        }
    }

    void scale(double s) {
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] * s;
        }
    }

    void translate(const Vector& t) {
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] + t;
        }
    }

	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
    Vector rho;
    bool mirror;
    BoundingBox bounding_box;
    BVHNode* root;
};

const Vector camera_center(0, 0, 55);
const Vector light_source(-10, 20, 40);
std::default_random_engine engine(1234);
std::uniform_real_distribution<double> unif_dist(0.0, 1.0);
const int width = 512;
const int height = 512;
const long long source_intensity = 3*1e8;
const double alpha = 60. * M_PI/180.;

class Sphere: public Object {
public:

    Sphere(const Vector &C, const double &R, const Vector &rho, bool mirror = false) : center(C), radius(R), ::Object(rho, mirror) {};

    bool intersect(const Ray &r, Intersection &I) const {
        double b = dot(r.direction, r.origin - center);
        double c = (r.origin - center).norm2() - radius * radius;
        double delta = b * b - c;
        if (delta >= 0) {
            double t1 = (dot(r.direction, center - r.origin) - sqrt(delta));
            double t2 = (dot(r.direction, center - r.origin) + sqrt(delta));
            if (t2 < 0) return false;
            I.t = (t1 > 0) ? t1 : t2;
            I.P = r.origin + I.t * r.direction;
            I.N = I.P - center;
            I.N.normalize();
            I.rho = rho;
            I.isMirror = mirror;
            return true;
        }
        return false;
    }

    Vector center;
    double radius;
};

class Scene {
public:
    std::vector<const Object*> objects;
    void addObject(const Object *s) {
        objects.push_back(s);
    }

    bool intersect(const Ray &r, Intersection &closest_intersection) {
        double min_t = std::numeric_limits<double>::max();
        bool intersection_found = false;

        for (const auto &object : objects) {
            Intersection I;
            if (object->intersect(r, I) && I.t < min_t) {
                min_t = I.t;
                closest_intersection = I;
                intersection_found = true;
            }
        }
        return intersection_found;
    }

    int visibility(const Vector &point, const Vector &light_source) {
        Vector light_direction = light_source - point;
        double light_distance = light_direction.norm();
        light_direction.normalize();

        Ray shadow_ray(point + 1e-6 * light_direction, light_direction);
        Intersection shadow_intersection;
        bool intersection_occurred = intersect(shadow_ray, shadow_intersection);

        if (intersection_occurred && shadow_intersection.t < light_distance) {
            return 0;
        } else {
            return 1;
        }
    }

    Vector getColor(const Ray &ray, int ray_depth, int max_path_length) {
        if (ray_depth == 0) return Vector(0., 0., 0.);

        Intersection closest_intersection;
        bool intersection_occurred = intersect(ray, closest_intersection);

        if (!intersection_occurred) {
            return Vector(0., 0., 0.);
        } else {
            if (closest_intersection.isMirror) {
                Vector reflected_dir = ray.direction - 2 * dot(ray.direction, closest_intersection.N) * closest_intersection.N;
                Ray reflected_ray(closest_intersection.P + 1e-4 * reflected_dir, reflected_dir);
                return getColor(reflected_ray, ray_depth - 1, max_path_length);
            } else {
                int visibility_term = visibility(closest_intersection.P, light_source);
                Vector color = (source_intensity / (4.0 * M_PI * (light_source - closest_intersection.P).norm2())) *
                            visibility_term * closest_intersection.rho * std::max(0.0, dot(closest_intersection.N, light_source - closest_intersection.P)) / M_PI;

                double r1 = unif_dist(engine);
                double r2 = unif_dist(engine);
                double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
                double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
                double z = sqrt(r2);

                Vector T1, T2;
                if (std::abs(closest_intersection.N[0]) < (std::abs(closest_intersection.N[1])) && std::abs(closest_intersection.N[0]) < (std::abs(closest_intersection.N[2]))) {
                    T1 = Vector(0, closest_intersection.N[2], -closest_intersection.N[1]);
                    T2 = Vector(1, 0, 0);
                } else if (std::abs(closest_intersection.N[1]) < (std::abs(closest_intersection.N[0])) && std::abs(closest_intersection.N[1]) < (std::abs(closest_intersection.N[2]))) {
                    T1 = Vector(closest_intersection.N[2], 0, -closest_intersection.N[0]);
                    T2 = Vector(0, 1, 0);
                } else {
                    T1 = Vector(closest_intersection.N[1], -closest_intersection.N[0], 0);
                    T2 = Vector(0, 0, 1);
                }

                Vector random_vec = x * T1 + y * T2 + z * closest_intersection.N;
                random_vec.normalize();

                color += closest_intersection.rho * getColor(Ray(closest_intersection.P + 1e-4 * random_vec, random_vec), ray_depth - 1, max_path_length);
                return color;
            }
        }
    }

};


int main() {
    std::vector<unsigned char> image(width * height * 3, 0);

    Scene scene;
    // scene.addObject(new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), true));
    scene.addObject(new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0)));
    scene.addObject(new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1)));
    scene.addObject(new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0)));
    scene.addObject(new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1)));
    scene.addObject(new Sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0)));
    scene.addObject(new Sphere(Vector(1000, 0, 0), 940, Vector(0, 1, 1)));

    int max_path_length = 5;
    int rays_per_pixel = 64;

    TriangleMesh mesh(Vector(1., 1., 1.), false);
    mesh.readOBJ("cat.obj");
    mesh.scale(0.6);
    mesh.translate(Vector(0, -10, 0));
    mesh.root = mesh.build_bvh(0, mesh.indices.size(), 4);
    scene.addObject(&mesh);

    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for schedule(dynamic, 1) 
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            Vector accumulated_color(0, 0, 0);
            for (int ray = 0; ray < rays_per_pixel; ray++) {
                double random_x_offset = dis(gen);
                double random_y_offset = dis(gen);
                Vector ray_dir;
                ray_dir[0] = (j + random_y_offset) - width / 2.0 + 0.5;
                ray_dir[1] = -(i + random_x_offset) + height / 2.0 + 0.5;
                ray_dir[2] = -width / (2.0 * std::tan(alpha / 2.0));

                ray_dir.normalize();
                Ray r(camera_center, ray_dir);

                accumulated_color += scene.getColor(r, max_path_length, max_path_length);
            }

            Vector color = accumulated_color / rays_per_pixel;
            image[(i * width + j) * 3 + 0] = std::min(255., std::pow(color[0], 1 / 2.2));
            image[(i * width + j) * 3 + 1] = std::min(255., std::pow(color[1], 1 / 2.2));
            image[(i * width + j) * 3 + 2] = std::min(255., std::pow(color[2], 1 / 2.2));
        }
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    double elapsedSeconds = elapsed.count();

    std::cout << "Elapsed time: " << elapsedSeconds << " seconds " << std::endl;

    stbi_write_png("image.png", width, height, 3, &image[0], 0);
    return 0;
}