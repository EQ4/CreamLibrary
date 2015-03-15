
#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <limits>

#ifndef DEF_HOA_DEFS_LIGHT
#define DEF_HOA_DEFS_LIGHT

#define HOA_PI  3.14159265358979323846264338327950288
#define HOA_2PI 6.283185307179586476925286766559005
#define HOA_PI2 1.57079632679489661923132169163975144
#define HOA_PI4 0.785398163397448309615660845819875721
#define HOA_1DEG (HOA_PI / 180.)

#include "../c.library.h"

using namespace std;

#if (__cplusplus <= 199711L)
#define noexcept
#endif

//#define HOA_EPSILON __FLT_EPSILON__ * 100.

namespace voro
{
    typedef unsigned long ulong;
    
    struct Triangle;
    
    struct Point
    {
        double x;
        double y;
        double z;
        ulong index;
        vector<Point> neightbours;
        vector<Point> bounds;
        
        Point() noexcept :
        x(0), y(0.), z(0.)
        {
            
        }
        
        ~Point() noexcept
        {
            neightbours.clear();
        }
        
        Point(double _x, double _y, double _z, ulong index = 0) noexcept :
        x(_x), y(_y), z(_z), index(0)
        {
            ;
        }
        
        Point(Point const& other) noexcept :
        x(other.x), y(other.y), z(other.z), index(other.index)
        {
            ;
        }
        
        double length(Point const& other) const noexcept
        {
            return sqrt((other.x - x) * (other.x - x) + (other.y - y) * (other.y - y) + (other.z - z) * (other.z - z));
        }
        
        double length() const noexcept
        {
            return sqrt((x + x) * (x + x) + (y + y) * (y + y) + (z + z) * (z + z));
        }
        
        double lenght2() const noexcept
        {
            const double l = length();
            return l * l;
        }
        
        Point operator-(Point const& other) const noexcept
        {
            return Point(x - other.x, y - other.y, z - other.z);
        }
        
        Point operator+(Point const& other) const noexcept
        {
            return Point(x + other.x, y + other.y, z + other.z);
        }
        
        Point operator*(Point const& other) const noexcept
        {
            return Point(x * other.x, y * other.y, z * other.z);
        }
        
        Point operator*(double val) const noexcept
        {
            return Point(x * val, y * val, z * val);
        }
        
        Point& operator*=(double val) noexcept
        {
            x *= val; y *= val; z *= val;
            return *this;
        }
        
        Point operator/(double val) const noexcept
        {
            return Point(x / val, y / val, z / val);
        }
        
        bool operator==(Point const& other) const noexcept
        {
            return other.x == x && y == other.y && other.z == z;
        }
        
        bool operator!=(Point const& other) const noexcept
        {
            return other.x != x || y != other.y || other.z != z;
        }
        
        Point cross(Point const& other) const noexcept
        {
            return Point(other.y * z - other.z * y, other.z * x - other.x * z, other.x * y - other.y * x);
        }
        
        double dot(Point const& other) const noexcept
        {
            return x * other.x + y * other.y + z * other.z;
        }
        
        double dot() const noexcept
        {
            return x * x + y * y + z * z;
        }
        
        static Point fromPolar(const double r, const double a, const double t) noexcept
        {
            return Point(r * cos(a + HOA_PI2) * cos(t), r * sin(a + HOA_PI2) * cos(t), r * sin(t));
        }
        
        void normalize() noexcept
        {
            const double l = length();
            if(l)
            {
                const double f = (2. / length());
                x *= f; y *= f; z *= f;
            }
            else
            {
                x = 0; y = 0; z = 0;
            }
        }
        
        Point normalized() const noexcept
        {
            Point t = *this;
            t.normalize();
            return t;
        }
        
        double radius()  const noexcept
        {
            return sqrt(x*x + y*y + z*z);
        }
        
        double azimuth() const noexcept
        {
            if (x == 0 && y == 0)
                return 0;
            return atan2(y, x) - HOA_PI2;
        }
        
        double elevation()  const noexcept
        {
            if(z == 0)
                return 0;
            return asin(z / sqrt(x*x + y*y + z*z));
        }
        
        double greatCircleDistance(Point const& other)  const noexcept
        {
            const double az1 = azimuth();
            const double az2 = other.azimuth();
            const double el1 = elevation();
            const double el2 = other.elevation();
            const double a = sin((az2 - az1) * 0.5);
            const double e = sin((el2 - el1) * 0.5);
            return 2. * asin(sqrt(a * a + cos(az1) * cos(az2) * e * e));
        }
        
        void addNeighbour(Point const& p)
        {
            if(find(neightbours.begin(), neightbours.end(), p) == neightbours.end())
            {
                neightbours.push_back(p);
            }
        }
        
        void addBound(Point const& p)
        {
            bounds.push_back(p);
        }
        
        void rotateZ(const double _z) noexcept
        {
            const double cosAngle = cos(_z);
            const double sinAngle = sin(_z);
            const double rx = x * cosAngle - y * sinAngle;
            y = x * sinAngle + y * cosAngle;
            x = rx;
        }
        
        void rotateY(const double _y) noexcept
        {
            const double cosAngle = cos(_y);
            const double sinAngle = sin(_y);
            const double rx = x * cosAngle - z * sinAngle;
            z = x * sinAngle + z * cosAngle;
            x = rx;
        }
        
        void rotateX(const double _x) noexcept
        {
            const double cosAngle = cos(_x);
            const double sinAngle = sin(_x);
            const double ry = y * cosAngle - z * sinAngle;
            z = y * sinAngle + z * cosAngle;
            y = ry;
        }
        
        void clean()
        {
            for(ulong i = 1; i < bounds.size(); )
            {
                if(bounds[i].length(bounds[i-1]) < 10e-6)
                {
                    bounds.erase(bounds.begin()+i);
                }
                else
                {
                    i++;
                }
            }
            if(bounds.size() > 1)
            {
                if(bounds[0].length(bounds[bounds.size()-1]) < 10e-6)
                {
                    bounds.pop_back();
                }
            }
        }
        
        void computeView(const bool top = true)
        {
            bool valid = false;
            for(ulong i = 0; i < bounds.size(); i++)
            {
                if(bounds[i].z > 0.)
                {
                    valid = true;
                }
            }
            if(!valid || bounds.size() < 3)
            {
                bounds.clear();
            }
            else
            {
                ulong size = bounds.size();
                for(ulong i = 0; i < size;)
                {
                    const ulong p = i ? i-1 : size-1;
                    const ulong n = (i == size-1) ? 0 : i+1;
                    if(bounds[i].z < 0. && bounds[p].z >= 0. && bounds[n].z >= 0.)
                    {
                        const double dist = bounds[p].z / (bounds[p].z - bounds[i].z);
                        Point temp1 = (bounds[i] - bounds[p]) * dist + bounds[p];
                        temp1.z = 0.;
                        temp1.normalize();
                        
                        bounds[i] = (bounds[i] - bounds[n]) * dist + bounds[n];
                        bounds[i].z = 0.;
                        bounds[i].normalize();
                        bounds.insert(bounds.begin()+i, temp1);
                        size++;
                        i += 3;
                    }
                    else if(bounds[i].z < 0. && bounds[p].z >= 0.)
                    {
                        const double dist = bounds[p].z / (bounds[p].z - bounds[i].z);
                        Point temp = (bounds[i] - bounds[p]) * dist + bounds[p];
                        temp.z = 0.;
                        temp.normalize();
                        bounds.insert(bounds.begin()+i, temp);
                        size++;
                        i += 2;
                    }
                    else if(bounds[i].z < 0. && bounds[n].z >= 0.)
                    {
                        const double dist = bounds[n].z / (bounds[n].z - bounds[i].z);
                        Point temp = (bounds[i] - bounds[n]) * dist + bounds[n];
                        temp.z = 0.;
                        temp.normalize();
                        bounds.insert(bounds.begin()+n, temp);
                        size++;
                        i += 2;
                    }
                    else
                    {
                        i++;
                    }
                }
                size = bounds.size();
                for(ulong i = 0; i < size;)
                {
                    const ulong p = i ? i-1 : size-1;
                    const ulong n = (i == size-1) ? 0 : i+1;
                    if(bounds[i].z <= 0. && bounds[p].z <= 0. && bounds[n].z <= 0.)
                    {
                        bounds.erase(bounds.begin()+i);
                        size--;
                    }
                    else
                    {
                        i++;
                    }
                }
            }
        }
        
        void sort()
        {
            const double el = HOA_PI2 - elevation();
            const double az = azimuth();
            for(ulong i = 0; i < bounds.size(); i++)
            {
                bounds[i].rotateZ(-az);
                bounds[i].rotateX(el);
            }
            std::sort(bounds.begin(), bounds.end(), compareAzimuth);
            for(ulong i = 0; i < bounds.size(); i++)
            {
                bounds[i].rotateX(-el);
                bounds[i].rotateZ(az);
            }
        }
        
        static bool compareAzimuth(Point const& p1, Point const& p2) noexcept
        {
            return p1.azimuth() < p2.azimuth();
        }
        
        static bool compareElevation(Point const& p1, Point const& p2) noexcept
        {
            return p1.elevation() < p2.elevation();
        }
    };
    
    struct Triangle
    {
        Point   a;
        Point   b;
        Point   c;
        Point   p;
        double  r;
        double  d;
        
        Triangle(Point const& _a, Point const& _b, Point const& _c) :
        a(_a), b(_b), c(_c)
        {
            const Point ac = (c - a);
            const Point ab = (b - a);
            const Point t = ab.cross(ac);
            const double _d = (2. * t.lenght2());
            if(_d > 1e-6)
            {
                p = (((t.cross(ab) * ac.lenght2()) + (ac.cross(t) * ab.lenght2())) / _d + a);
                if(p.length(Point(0., 0., 0)) > 1e-6)
                {
                    p.normalize();
                    r = p.length(a);
                    d = p.greatCircleDistance(a);
                }
                else
                {
                    r = 0.;
                    d = 0.;
                }
            }
            else
            {
                r = 0.;
            }
        }
        
        ~Triangle()
        {
            ;
        }
    };
    
    class Voronoi
    {
    private:
        vector<Point>       m_points;
        vector<Triangle>    m_triangles;
    public:
        
        Voronoi() noexcept
        {
            ;
        }
        
        ~Voronoi() noexcept
        {
            clear();
        }
        
        void add(Point const& p)
        {
            m_points.push_back(p.normalized());
            m_points[m_points.size()-1].index = m_points.size();
        }
        
        void clear()
        {
            m_points.clear();
            m_triangles.clear();
        }

        vector<Point> const& getPoints() const noexcept
        {
            return m_points;
        }
        
        vector<Point>& getPoints() noexcept
        {
            return m_points;
        }
        
        void compute()
        {
            m_triangles.clear();
            std::sort(m_points.begin(), m_points.end(), Point::compareElevation);
            
            for(ulong i = 0; i < m_points.size() - 2; i++)
            {
                for(ulong j = i+1; j < m_points.size() - 1; j++)
                {
                    for(ulong k = j+1; k < m_points.size(); k++)
                    {
                        Triangle t(m_points[i], m_points[j], m_points[k]);
                        if(t.r != 0.)
                        {
                            bool valid = true;
                            for(ulong l = 0; l < m_points.size(); l++)
                            {
                                if(l != i && l != j && l != k)
                                {
                                    if(t.p.length(m_points[l]) < t.r - 1e-6)
                                    {
                                        valid = false;
                                    }
                                }
                            }
                            if(valid)
                            {
                                m_triangles.push_back(t);
                                m_points[i].addNeighbour(m_points[j]);
                                m_points[i].addNeighbour(m_points[k]);
                                m_points[i].addBound(t.p);
                                m_points[j].addNeighbour(m_points[i]);
                                m_points[j].addNeighbour(m_points[k]);
                                m_points[j].addBound(t.p);
                                m_points[k].addNeighbour(m_points[i]);
                                m_points[k].addNeighbour(m_points[j]);
                                m_points[k].addBound(t.p);
                            }
                        }
                    }
                }
            }
            for(ulong i = 0; i < m_points.size(); i++)
            {
                m_points[i].sort();
                m_points[i].clean();
                m_points[i].computeView();
                m_points[i].clean();
            }
        }
    };
}

#endif