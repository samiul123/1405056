#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include "bits/stdc++.h"
#include "bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))
#define epsilon (1.0e-6)
#define deg_to_rad (pi / 180)
#define INSIDE_RANGE 1
#define OUTSIDE_NEAR 2
#define OUTSIDE_FAR 3

class Vector
{
public:
	double x, y, z;

	// constructs a vector with given components
	Vector(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	// keeps the direction same. recalculates the vector to be unit.
	void normalize()
	{
		double r = sqrt(x*x + y * y + z * z);
		x = x / r;
		y = y / r;
		z = z / r;
	}

	// add two vectors
	Vector operator+(const Vector& v)
	{
		Vector v1(x + v.x, y + v.y, z + v.z);
		return v1;
	}

	// subtract one vector from another
	Vector operator-(const Vector& v)
	{
		Vector v1(x - v.x, y - v.y, z - v.z);
		return v1;
	}

	// scale a vector with a given coefficient
	Vector operator* (double m)
	{
		Vector v(x*m, y*m, z*m);
		return v;
	}

	// get the dot product of two vectors
	static double dot(Vector a, Vector b)
	{
		return a.x*b.x + a.y*b.y + a.z*b.z;
	}

	// get the cross product of two vectors
	static Vector cross(Vector a, Vector b)
	{
		Vector v(a.y*b.z - a.z*b.y, b.x*a.z - b.z*a.x, a.x*b.y - a.y*b.x);
		return v;
	}

	// print a vector. only for testing purposes.
	void print()
	{
		cout << "Vector" << endl;
		cout << x << " " << y << " " << z << endl;
	}
};

class homogeneous_point
{
public:
	double x, y, z, w;

	// set the three coordinates, set w to 1
	homogeneous_point(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = 1;
	}

	/*
	default constructor. does nothing. allows declarations like below:
		matrix m;
	therefore, usage is dangerous
	*/
	homogeneous_point() {
	}

	// constructs a homogeneous point with given coordinates. forces w to be 1.0
	// if w is zero, raises error
	homogeneous_point(double x, double y, double z, double w)
	{
		assert(w != 0);
		this->x = x / w;
		this->y = y / w;
		this->z = z / w;
		this->w = 1;
	}

	void operator=(const homogeneous_point point) {
		this->x = point.x;
		this->y = point.y;
		this->z = point.z;
		this->w = point.w;
	}
	// adds two points. returns a point forcing w to be 1.0
	homogeneous_point operator+ (const homogeneous_point& point)
	{
		double x = this->x + point.x;
		double y = this->y + point.y;
		double z = this->z + point.z;
		double w = this->w + point.w;
		homogeneous_point p(x, y, z, w);
		return p;
	}

	// subtracts one point from another. returns a point forcing w to be 1.0
	homogeneous_point operator- (const homogeneous_point& point)
	{
		double x = this->x - point.x;
		double y = this->y - point.y;
		double z = this->z - point.z;
		double w = this->w - point.w;
		homogeneous_point p(x, y, z, w);
		return p;
	}

	bool operator==(const homogeneous_point point) {
		return this->x == point.x && this->y == point.y && this->z == point.z;
	}

	// Print the coordinates of a point. exists for testing purpose.
	void print()
	{
		cout << "Point: " << endl;
		cout << x << " " << y << " " << z << " " << w << endl;
	}

	//point + vector
	homogeneous_point operator+ (const Vector& vector)
	{
		double x = this->x + vector.x;
		double y = this->y + vector.y;
		double z = this->z + vector.z;
		double w = this->w;
		homogeneous_point p(x, y, z, w);
		return p;
	}
};




/*
The matrices are forced to be 4x4. This is because in this assignment, we will deal with points in triangles.
Maximum # of points that we will deal with at once is 3. And all the standard matrices are 4x4 (i.e. scale, translation, rotation etc.)
*/
class matrix
{
public:
	double values[4][4];
	int num_rows, num_cols;

	// only set the number of rows and cols
	matrix(int rows, int cols)
	{
		assert(rows <= 4 && cols <= 4);
		num_rows = rows;
		num_cols = cols;
	}

	// prepare an nxn square matrix
	matrix(int n)
	{
		assert(n <= 4);
		num_rows = num_cols = n;
	}

	// prepare and return an identity matrix of size nxn
	static matrix make_identity(int n)
	{
		assert(n <= 4);
		matrix m(n);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)
					m.values[i][j] = 1;
				else
					m.values[i][j] = 0;
			}
		}
		return m;
	}

	// print the matrix. exists for testing purposes
	void print()
	{
		cout << "Matrix:" << endl;
		for (int i = 0; i < num_rows; i++)
		{
			for (int j = 0; j < num_cols; j++)
			{
				cout << values[i][j] << "\t";
			}
			cout << endl;
		}
	}

	// add the two matrices. Raise error if dimension mismatches
	matrix operator+ (const matrix& m)
	{
		assert(this->num_rows == m.num_rows);
		assert(this->num_cols == m.num_cols);

		matrix m1(num_rows, num_cols);
		for (int i = 0; i < num_rows; i++)
		{
			for (int j = 0; j < num_cols; j++)
			{
				m1.values[i][j] = values[i][j] + m.values[i][j];
			}
		}
		return m1;
	}

	// subtract a matrix from another. raise error if dimension mismatches
	matrix operator- (const matrix& m)
	{
		assert(this->num_rows == m.num_rows);
		assert(this->num_cols == m.num_cols);

		matrix m1(num_rows, num_cols);
		for (int i = 0; i < num_rows; i++)
		{
			for (int j = 0; j < num_cols; j++)
			{
				m1.values[i][j] = values[i][j] - m.values[i][j];
			}
		}
		return m1;
	}

	// multiply two matrices. allows statements like m1 = m2 * m3; raises error is dimension mismatches
	matrix operator* (const matrix& m)
	{
		assert(this->num_cols == m.num_rows);
		matrix m1(this->num_rows, m.num_cols);

		for (int i = 0; i < m1.num_rows; i++) {
			for (int j = 0; j < m1.num_cols; j++) {
				double val = 0;
				for (int k = 0; k < this->num_cols; k++) {
					val += this->values[i][k] * m.values[k][j];
				}
				m1.values[i][j] = val;
			}
		}
		return m1;
	}

	// multiply a matrix with a constant
	matrix operator* (double m)
	{
		matrix m1(this->num_rows, this->num_cols);
		for (int i = 0; i < num_rows; i++) {
			for (int j = 0; j < num_cols; j++) {
				m1.values[i][j] = m * this->values[i][j];
			}
		}
		return m1;
	}

	// multiply a 4x4 matrix with a homogeneous point and return the resulting point.
	// usage: homogeneous_point p = m * p1;
	// here, m is a 4x4 matrix, intended to be the transformation matrix
	// p1 is the point on which the transformation is being made
	// p is the resulting homogeneous point
	homogeneous_point operator* (const homogeneous_point& p)
	{
		assert(this->num_rows == this->num_cols && this->num_rows == 4);

		matrix m(4, 1);
		m.values[0][0] = p.x;
		m.values[1][0] = p.y;
		m.values[2][0] = p.z;
		m.values[3][0] = p.w;

		matrix m1 = (*this)*m;
		homogeneous_point p1(m1.values[0][0], m1.values[1][0], m1.values[2][0], m1.values[3][0]);
		return p1;
	}

	// return the transpose of a matrix
	matrix transpose()
	{
		matrix m(num_cols, num_rows);
		for (int i = 0; i < num_rows; i++) {
			for (int j = 0; j < num_cols; j++) {
				m.values[j][i] = values[i][j];
			}
		}
		return m;
	}

};

/*
A simple class to hold the color components, r, g, b of a certain shade.
*/
class color {
public:
	double r, g, b;
	color(double r, double g, double b) {
		this->r = r;
		this->g = g;
		this->b = b;
	}
	color() {
	}
	void operator=(const color& c) {
		this->r = c.r;
		this->g = c.g;
		this->b = c.b;
	}
	void print() {
		cout << "color" << endl;
		cout << r << " " << g << " " << b << endl;
	}
};

class Triangle {
public:
	homogeneous_point p1, p2, p3;
	color c;
	Triangle(homogeneous_point p1, homogeneous_point p2, homogeneous_point p3, color c) {
		this->p1 = p1;
		this->p2 = p2;
		this->p3 = p3;
		this->c = c;
	}
	Triangle(homogeneous_point p1, homogeneous_point p2, homogeneous_point p3) {
		this->p1 = p1;
		this->p2 = p2;
		this->p3 = p3;
	}
	void operator=(const Triangle& triangle) {
		this->p1 = triangle.p1;
		this->p2 = triangle.p2;
		this->p3 = triangle.p3;
	}
	void print() {
		cout << "Triangle:" << endl;
		p1.print(); 
		p2.print();
		p3.print();
		c.print();
	}
};




double eye_x, eye_y, eye_z;
double look_x, look_y, look_z;
double up_x, up_y, up_z;
double fov_x, fov_y, aspectRatio, near, far;
color backgroud;
int screen_x, screen_y;
vector<homogeneous_point> points;
vector<color> point_color;
vector<Triangle> triangles;
vector<Triangle> stage3_triangles;
map<pair<double, double>, int> map_x;
map<pair<double, double>, int> map_y;

double str_to_double(string s) {
	double u = 0;
	stringstream geek(s);
	geek >> u;
	return u;
}

bool compare(homogeneous_point p1, homogeneous_point p2) {
	return p1.y < p2.y;
}

void map_x_y() {
	double diff_x = 2.0 / screen_x;
	double diff_y = 2.0 / screen_y;
	int i = 0;
	double sum_x = -1;
	double prev_x = 0;
	while (sum_x < 1) {
		prev_x = sum_x;
		sum_x += diff_x;
		map_x[make_pair(prev_x, sum_x)] = i;
		i++;
	}
	i = 0;
	double sum_y = -1;
	double prev_y = 0;
	while (sum_y < 1) {
		prev_y = sum_y;
		sum_y += diff_y;
		map_y[make_pair(prev_y, sum_y)] = i;
		i++;
	}
}

homogeneous_point find_intersection_point(homogeneous_point inside_point, homogeneous_point outside_point
	, int range_indicator) {
	homogeneous_point diff_point = outside_point - inside_point;
	Vector diff_vector(diff_point.x, diff_point.y, diff_point.z);
	double t;
	if (range_indicator == OUTSIDE_NEAR)t = (-near - inside_point.z) / diff_point.z;
	else if (range_indicator == OUTSIDE_FAR)t = (-far - inside_point.z) / diff_point.z;
	homogeneous_point intersected_point = inside_point + diff_vector * t;
	cout << "Intersection point:" << endl;
	intersected_point.print();
	return intersected_point;
}

void print_triangle(vector<Triangle> triangles) {
	for (auto t : triangles) {
		t.print();
	}
}

int is_in_range(double z, int i) {
	/*if (z <= -near && z >= -far) return INSIDE_RANGE;
	else if (z > -near) return OUTSIDE_NEAR;
	else if (z < -far)return OUTSIDE_FAR;
	return 0;*/
	if (i == 0) {
		if (z <= -near) {
			cout << "i: " << i << " " << z << " inside" << endl;
			return INSIDE_RANGE;
		}
		else if (z > -near) {
			cout << "i: " << i << " " << z << " outside near" << endl;
			return OUTSIDE_NEAR;
		}
	}
	else {
		if (z >= -far) {
			cout << "i: " << i << " " << z << " inside" << endl;
			return INSIDE_RANGE;
		}
		else if (z < -far) {
			cout << "i: " << i << " " << z << " outside far" << endl;
			return OUTSIDE_FAR;
		}
	}
	return 0;
}



void scan_convert() {
	ifstream stage3;
	stage3.open("stage3.txt");

	color** pixels = new color*[screen_x];
	double** zs = new double*[screen_x];
	for (int i = 0; i < screen_x; i++) {
		pixels[i] = new color[screen_y];
		for (int j = 0; j < screen_y; j++) {
			pixels[i][j] = backgroud;
		}
		zs[i] = new double[screen_y];
		for (int j = 0; j < screen_y; j++) {
			zs[i][j] = +20; // a very large value intended as +INFINITY
		}
	}
	
	map_x_y();
	for (auto t: stage3_triangles) {
		homogeneous_point point_arr[] = { t.p1, t.p2, t.p3 };
		cout << "before sorting" << endl;
		for (int i = 0; i < 3; i++) {
			point_arr[i].print();
		}
		cout << endl;
		sort(point_arr, point_arr + 3, compare);
		cout << "after sorting" << endl;
		for (int i = 0; i < 3; i++) {
			point_arr[i].print();
		}
		cout << endl;
	}
	// perform scan conversion, populate the 2D array pixels
	// the array zs is the z-buffer.


	// the following code generates a bmp image. do not change this.
	bitmap_image image(screen_x, screen_y);
	for (int x = 0; x < screen_x; x++) {
		for (int y = 0; y < screen_y; y++) {
			image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
		}
	}
	image.save_image("out.bmp");

	// free the dynamically allocated memory

}




void stage3()
{
	if (near == far) return;
	ifstream stage2;
	ofstream stage3;
	stage2.open("stage2.txt");
	stage3.open("stage3.txt");
	stage3 << std::fixed;
	stage3 << std::setprecision(7);

	// process input from stage2 and write to stage3
	fov_x = fov_y * aspectRatio;
	double t = near * tan((fov_y / 2)*deg_to_rad);
	double r = near * tan((fov_x / 2)*deg_to_rad);
	matrix projection_matrix = matrix::make_identity(4);
	projection_matrix.values[0][0] = near / r;
	projection_matrix.values[1][1] = near / t;
	projection_matrix.values[2][2] = -(far + near) / (far - near);
	projection_matrix.values[3][2] = -1;
	projection_matrix.values[2][3] = -(2*far*near)/(far - near);
	projection_matrix.values[3][3] = 0;

	string line;
	int p = 0;
	int q = 0;
	vector<homogeneous_point> points;
	while (getline(stage2, line)) {
		istringstream iss(line);
		vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
		if (!results.size()) continue;
		homogeneous_point point(str_to_double(results[0]), str_to_double(results[1]), str_to_double(results[2]));
		/*if (point.z > -near)point.z = -near;
		else if(point.z < -far)point.z = -far;*/
/*
		
		homogeneous_point transformed_point = projection_matrix * point;
		stage3 << transformed_point.x << " " << transformed_point.y << " " << transformed_point.z << endl;*/
		points.push_back(point);
		p++;
		if (p == 3) {
			//stage3 << endl;
			//cout << points.size() << endl;
			int clip = 0;
			color c = triangles[q].c;
			q++;
			for (int i = 0; i <= 1; i++) {
				homogeneous_point s = points[(points.size() - 1)%3];
				cout << "s: ";
				s.print();
				homogeneous_point p = points[0];
				cout << "p: ";
				p.print();
				int j = 0;
				vector<homogeneous_point> outList;
				cout << "j: " << j << endl;
				cout << "points size: " << points.size() << endl;
				for (auto p:points) {
					p.print();
				}
				while (j < points.size()) {
					if (is_in_range(s.z, i) == INSIDE_RANGE && is_in_range(p.z, i) == INSIDE_RANGE) outList.push_back(p);
					else if (is_in_range(s.z, i) == INSIDE_RANGE) {
						clip = 1;
						if (is_in_range(p.z, i) == OUTSIDE_NEAR) {
							homogeneous_point intersection_point = find_intersection_point(s, p, OUTSIDE_NEAR);
							outList.push_back(intersection_point);
						}
						else if (is_in_range(p.z, i) == OUTSIDE_FAR) {
							homogeneous_point intersection_point = find_intersection_point(s, p, OUTSIDE_FAR);
							outList.push_back(intersection_point);
						}
					}
					else if (is_in_range(p.z, i) == INSIDE_RANGE) {
						clip = 1;
						if (is_in_range(s.z, i) == OUTSIDE_NEAR) {
							homogeneous_point intersection_point = find_intersection_point(p, s, OUTSIDE_NEAR);
							outList.push_back(intersection_point);
							outList.push_back(p);
						}
						else if (is_in_range(s.z, i) == OUTSIDE_FAR) {
							homogeneous_point intersection_point = find_intersection_point(p, s, OUTSIDE_FAR);
							outList.push_back(intersection_point);
							outList.push_back(p);
						}
						
					}
					s = points[j];
					p = points[(j + 1)%3];
					j++;
				}
				
				points.clear();
				for (auto t:outList) {
					points.push_back(t);
				}
			}
			cout << "after clip: " << endl;
			for (auto t:points) {
				t.print();
			}
			//cout << points.size() << endl;
			if (clip) {
				homogeneous_point clipped_point_1 = points[0];
				int m = 1;
				while (m < (points.size() - 1)) {
					homogeneous_point clipped_point_2 = points[m];
					homogeneous_point clipped_point_3 = points[m + 1];
					Triangle triangle(clipped_point_1, clipped_point_2, clipped_point_3, c);
					stage3_triangles.push_back(triangle);
					m++;
				}
			}
			else {
				Triangle triangle(points[0], points[1], points[2], c);
				stage3_triangles.push_back(triangle);
			}
			p = 0;
			points.clear();
		}
	}
	print_triangle(stage3_triangles);
	for (auto t:stage3_triangles) {
		t.p1 = projection_matrix * t.p1;
		t.p2 = projection_matrix * t.p2;
		t.p3 = projection_matrix * t.p3;
		/*t.p1 = p1;
		t.p2 = p2;
		t.p3 = p3;*/
		stage3 << t.p1.x << " " << t.p1.y << " " << t.p1.z << endl;
		stage3 << t.p2.x << " " << t.p2.y << " " << t.p2.z << endl;
		stage3 << t.p3.x << " " << t.p3.y << " " << t.p3.z << endl << endl;
	}
	
	stage3.close();
	stage2.close();

}

void stage2()
{
	ifstream stage1;
	ofstream stage2;
	stage1.open("stage1.txt");
	stage2.open("stage2.txt");
	stage2 << std::fixed;
	stage2 << std::setprecision(7);

	// collect input from stage1 and process, write output to stage2
	Vector look(look_x, look_y, look_z);
	Vector eye(eye_x, eye_y, eye_z);
	Vector l = look - eye;
	l.normalize();
	Vector up(up_x, up_y, up_z);
	Vector r = Vector::cross(l, up);
	r.normalize();
	Vector u = Vector::cross(r, l);

	matrix R = matrix::make_identity(4);
	R.values[0][0] = r.x;
	R.values[0][1] = r.y;
	R.values[0][2] = r.z;
	R.values[1][0] = u.x;
	R.values[1][1] = u.y;
	R.values[1][2] = u.z;
	R.values[2][0] = -l.x;
	R.values[2][1] = -l.y;
	R.values[2][2] = -l.z;

	matrix T = matrix::make_identity(4);
	T.values[0][3] = -eye_x;
	T.values[1][3] = -eye_y;
	T.values[2][3] = -eye_z;

	matrix V = R * T;
	string line;
	int i = 0;
	int j = 0;
	vector<homogeneous_point> points;
	while (getline(stage1, line)){
		istringstream iss(line);
		vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
		if (!results.size()) continue;
		homogeneous_point point(str_to_double(results[0]), str_to_double(results[1]), str_to_double(results[2]));
		homogeneous_point transformed_point = V * point;
		points.push_back(transformed_point);
		stage2 << transformed_point.x << " "<<  transformed_point.y << " " << transformed_point.z << endl;
		i++;
		if (i == 3) {
			Triangle triangle(points[0], points[1], points[2]);
			triangles[j] = triangle;
			j++;
			stage2 << endl;
			i = 0;
		}
	}
	//print_triangle(triangles);
	stage1.close();
	stage2.close();

}

Vector rodrigues_formula(Vector x, Vector a, double angle) {
	//homogeneous_point p;
	Vector p = x * cos(angle * deg_to_rad) + a * (1 - cos(angle * deg_to_rad))*(Vector::dot(x, a)) + Vector::cross(a, x)*sin(angle*deg_to_rad);
	return p;
}

void stage1()
{
	ifstream scene;
	ofstream stage1;
	scene.open("scene.txt");
	stage1.open("stage1.txt");
	stage1 << std::fixed;
	stage1 << std::setprecision(7);

	string command;

	scene >> eye_x >> eye_y >> eye_z;
	scene >> look_x >> look_y >> look_z;
	scene >> up_x >> up_y >> up_z;
	scene >> fov_y >> aspectRatio >> near >> far;
	scene >> screen_x >> screen_y;
	scene >> backgroud.r >> backgroud.g >> backgroud.b;

	// take other commands as input from scene in a loop
	// process accordingly
	// write to stage1
	stack<matrix> s;
	//matrix i_matrix(4);
	s.push(matrix::make_identity(4));
	while (true) {
		scene >> command;
		if (command == "triangle") {
			homogeneous_point points[3];
			int r, g, b;
			scene >> points[0].x >> points[0].y >> points[0].z;
			scene >> points[1].x >> points[1].y >> points[1].z;
			scene >> points[2].x >> points[2].y >> points[2].z;
			scene >> r >> g >> b;
			color c(r, g, b);
			Triangle triangle(points[0], points[1], points[2], c);
			triangles.push_back(triangle);
			for (auto x: points) {
				x.w = 1;
				x = s.top() * x;
				stage1 << x.x << " " << x.y << " " << x.z << endl;
			}
			stage1 << endl << endl;
		}
		else if (command == "translate") {
			double t_x, t_y, t_z;
			scene >> t_x >> t_y >> t_z;
			matrix translation_matrix = matrix::make_identity(4);
			translation_matrix.values[0][3] = t_x;
			translation_matrix.values[1][3] = t_y;
			translation_matrix.values[2][3] = t_z;
			//s.push(translation_matrix*s.top());
			s.top() = s.top() * translation_matrix;
		}
		else if (command == "scale") {
			double s_x, s_y, s_z;
			scene >> s_x >> s_y >> s_z;
			matrix scale_matrix = matrix::make_identity(4);
			scale_matrix.values[0][0] = s_x;
			scale_matrix.values[1][1] = s_y;
			scale_matrix.values[2][2] = s_z;
			//s.push(scale_matrix*s.top());
			s.top() = s.top() * scale_matrix;
		}
		else if (command == "rotate") {
			double angle, a_x, a_y, a_z;
			scene >> angle >> a_x >> a_y >> a_z;
			Vector axis_a(a_x, a_y, a_z);
			axis_a.normalize();
			Vector i(1, 0, 0);
			Vector j(0, 1, 0);
			Vector k(0, 0, 1);
			Vector c1 = rodrigues_formula(i, axis_a, angle);
			Vector c2 = rodrigues_formula(j, axis_a, angle);
			Vector c3 = rodrigues_formula(k, axis_a, angle);
			matrix rotation_matrix = matrix::make_identity(4);
			rotation_matrix.values[0][0] = c1.x;
			rotation_matrix.values[1][0] = c1.y;
			rotation_matrix.values[2][0] = c1.z;
			rotation_matrix.values[0][1] = c2.x;
			rotation_matrix.values[1][1] = c2.y;
			rotation_matrix.values[2][1] = c2.z;
			rotation_matrix.values[0][2] = c3.x;
			rotation_matrix.values[1][2] = c3.y;
			rotation_matrix.values[2][2] = c3.z;
			s.top()  = s.top()*rotation_matrix;
			/*s.pop();
			s.push(t);*/
			//s.push(rotation_matrix*s.top());
		}
		else if (command == "push") {
			s.push(s.top());
		}
		else if (command == "pop") {
			s.pop();
		}
		else if (command == "end") {
			break;
		}
	}
	



	scene.close();
	stage1.close();

}

int main()
{
	cout << std::fixed;
	cout << std::setprecision(4);

	stage1();
	stage2();
	stage3();
	scan_convert();
	cin.get();
	return 0;
}
