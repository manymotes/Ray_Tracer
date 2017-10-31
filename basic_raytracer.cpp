#include <fstream>
#include <cmath>
#include <string>
#include <iostream>


struct Vec3 {
    double x,y,z;
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
    Vec3 operator + (const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
    Vec3 operator - (const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
    Vec3 operator * (double d) const { return Vec3(x*d, y*d, z*d); }
    Vec3 operator / (double d) const { return Vec3(x/d, y/d, z/d); }
    Vec3 normalize() const {
        double mg = sqrt(x*x + y*y + z*z);
        return Vec3(x/mg,y/mg,z/mg);
        
    
    }
    double findLength() const
    {
        double toreturn =sqrt(x*x + y*y + z*z);
        return toreturn;
    }
};


inline double dot(const Vec3& a, const Vec3& b) {
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}
Vec3 MultiplyScalar(Vec3 a, float b) { return Vec3(a.x * b, a.y * b, a.z * b); }

bool sign(int input) {
    if (input > 0.1f)
        return true;
    else
        return false;
}

static Vec3 CrossProduct(Vec3 v1, Vec3 v2)
{
    
    Vec3 toReturn = Vec3
    (
     v1.y * v2.z - v1.z * v2.y,
     v1.z * v2.x - v1.x * v2.z,
     v1.x * v2.y - v1.z * v2.x
     );
    
    return toReturn;
}

struct Ray {
    Vec3 o,d;
    Ray(const Vec3& o, const Vec3& d) : o(o), d(d) {}
};

//struct Sphere {
//    Vec3 c;
//    double r;
//    Sphere(const Vec3& c, double r) : c(c), r(r) {}
//    Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
//    bool intersect(const Ray& ray, double &t) const {
//        const Vec3 o = ray.o;
//        const Vec3 d = ray.d;
//        const Vec3 oc = c - o;
//        const double tca = dot(d.normalize(), oc);
//
//        if (tca < 0)
//        {
//            return false;
//        }
//
//        const double thc = r*r - oc.findLength() + tca * tca;
//
//        if (thc < 0)
//        {
//            return false;
//        }
//
//        t = tca - thc;
//
//
//
////        const double b = 2 * dot(oc, d);
////        const double c = dot(oc, oc) - r*r;
////        double disc = b*b - 4 * c;
//
//
////        if (disc < 1e-4) return false;
////        disc = sqrt(disc);
////        const double t0 = -b - disc;
////        const double t1 = -b + disc;
////        t = (t0 < t1) ? t0 : t1;
//
//
//
//        return true;
//    }
//};


struct Sphere {
    Vec3 c;
    double r;
    Sphere(const Vec3& c, double r) : c(c), r(r) {}
    Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
    
    
    bool intersect(const Ray& ray, double &t) const {
        const Vec3 o = ray.o;
        const Vec3 d = ray.d;
        const Vec3 oc = o - c;
        const double b = 2 * dot(oc, d);
        const double c = dot(oc, oc) - r*r;
        double disc = b*b - 4 * c;
        if (disc < 1e-4) return false;
        disc = sqrt(disc);
        const double t0 = -b - disc;
        const double t1 = -b + disc;
        
        const double tca = dot(d.normalize(), oc);
        const double thc = r*r - oc.findLength() + tca * tca;
       // t = tca - thc;
        t = (t0 < t1) ? t0 : t1;
        return true;
    }
};


struct Sphere2 {
    Vec3 c;
    double r;
    Sphere2(const Vec3& c, double r) : c(c), r(r) {}
    Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
   
    
    bool intersect(const Ray& ray, double &t) const {
        const Vec3 o = ray.o;
        const Vec3 d = ray.d;
        const Vec3 oc = o - c;
        const double b = 2 * dot(oc, d);
        const double c = dot(oc, oc) - r*r;
        double disc = b*b - 4 * c;
        if (disc < 1e-4) return false;
        disc = sqrt(disc);
        const double t0 = -b - disc;
        const double t1 = -b + disc;
        
        const double tca = dot(d.normalize(), oc);
        const double thc = r*r - oc.findLength() + tca * tca;
        // t = tca - thc;
        t = (t0 < t1) ? t0 : t1;
        return true;
    }
};
struct Triangle {
    Vec3 a;
    Vec3 b;
    Vec3 c;
    
    
    //float TriangleDistance;
    
    
    Triangle(const Vec3& a, const Vec3& b, const Vec3& c) : a(a), b(b), c(c)  {
    }
    //  Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
    
//    float DistanceToTriangle()
//    {
//        return TriangleDistance;
//    }
    Vec3 getNormal() const
    {
        Vec3 v1 = a-b;
        Vec3 v2 = c-b;
        Vec3 pn = CrossProduct(v2, v1);
        return pn.normalize();
    }
    bool intersect(const Ray& ray, double &t_out) const {
        
       
        Vec3 pn = getNormal();
        
        const Vec3 o = ray.o;
        const Vec3 d = ray.d;
        //Vec3 normal = CrossProduct(b-a,tA),(b-c)).normalize();
        double dis =  -dot(pn, b);
        
        double t = -((pn.x * o.x) + (pn.y * o.y) + (pn.z * o.z) + dis)/( (pn.x * d.x) + (pn.y * d.y) + (pn.z * d.z));
        
        Vec3 intersect = (o + MultiplyScalar(d, t));
        
        bool signA = sign(dot(CrossProduct((b - a),(intersect - a)),pn));
        bool signB = sign(dot(CrossProduct((c - b),(intersect - b)),pn));
        bool signC = sign(dot(CrossProduct((a - c),(intersect - c)),pn));
        
        
        //help with this t value
        //t_out = (t0 < t1) ? t0 : t1;
        t_out = t;
        return ((signA == signB) && (signA == signC) && t>0.1);
        
        
        
        
        
        
//        const Vec3 o = ray.o;
//        const Vec3 d = ray.d;
//        //        const Vec3 oc = o - c;
//        //        const double b = 2 * dot(oc, d);
//        //        const double c = dot(oc, oc) - r*r;
//        //        double disc = b*b - 4 * c;
//        Vec3 v1 = a-b;
//        Vec3 v2 = b-c;
//        Vec3 pn = CrossProduct(a, b);
//
//        ////!!!!! find dominant axis!!!!!
//        std::string dominatAxis;
//        if (abs(pn.x) > abs(pn.y))
//        {
//            if(abs(pn.x) >= abs(pn.z))
//            {
//                dominatAxis = "x";
//            }
//            else{
//                dominatAxis ="z";
//            }
//        }
//        else{
//            if (abs(pn.y) >= abs(pn.z))
//            {
//                dominatAxis = "y";
//            }
//            else
//            {
//                dominatAxis = "z";
//            }
//        }
//
//        //gets intersection of ray and plane
//
//        double dis = -dot(pn, b);
//        Vec3 intersect = -((pn.x * o.x) + (pn.y * o.y) + (pn.z * o.z) + dis)/( (pn.x * d.x) + (pn.y * d.y) + (pn.z * d.z));
//
//
//        Vecw projectedVertices[3];
//        Vec2 ri;
//
//        if (dominatAxis == "x")
//        {
//            projectedVertices[0] = Vec2(a.y, a.z);
//            projectedVertices[1] =Vec2(b.y, b.z);
//            projectedVertices[2] = Vec2(c.y, c.z);
//            ri = Vec2(intersect.y, intersect.z);
//        }
//        else if (dominatAxis == "y")
//        {
//            projectedVertices[0] = Vec2(a.x, a.z);
//            projectedVertices[1] = Vec2(b.x, b.z);
//            projectedVertices[2] = Vec2(c.x, c.z);
//            ri = Vec2(intersect.x, intersect.z);
//        }
//        else if(dominatAxis == "z")
//        {
//            projectedVertices[0] = Vec2(a.x, a.y);
//            projectedVertices[1] = Vec2(b.x, b.y);
//            projectedVertices[2] = Vec2(c.x, c.y);
//            ri = Vec2(intersect.x, intersect.y);
//        }
//
//
//
//        //     ????? how do you find d
//
//
//
//        //translate verticies by -ri
//        g0 = g0 - ri;
//        g1 = g1 - ri;
//        g2 = g2 - ri;
//        Vec2 trianlgeEdges[3] = {g0, g1, g2};
//
//        ri = ri - ri;
//        if (ri != Vec2(0,0))
//        {
//            throw noexcept;
//        }
//
//        //u and v prime are g0 g2 and g2
//        //singholde, look at the v coordinate of the first point. which is 3 in the example.
//        // sign holder = 1, becuase 3 is positive.
//        if v1 < 0 signHolder = -1;
//        else singholder = 1;
//
//        int nextsing;
//        int sumCrossing = 0;
//
//        for( int i  = 0, i < 2, i++)
//        {
//            if (trianlge[i+1] < 0)
//            {
//                nextsing = -1;
//            }
//            else {
//                nextsing = 1;
//            }
//
//            if (singholder != nextsing)
//            {
//                if(projectedVertices[i].x > 0 && projectedVertices[i+1].x > 0)
//                {
//                    numCrossing +=1;
//                }
//                else if(projectedVertices[i].x > 0 || projectedVertices[i+1].x > 0)
//                {
//                    float uCross = projectedVertices[i].x - projectedVertices[i].y (projectedVertices[i+1].x - projectedVertices[i].x)/(projectedVertices[i+1].y - projectedVertices[i].y);
//
//                    if (uCross > 0)
//                    {
//                        numCrossing += 1;
//                    }
//                }
//            }
//        }
//
//        if (numCrossing == 1)
//        {
//            return true;
//        }
//        return false;
//    }
}
};

void clamp255(Vec3& col) {
    
 
    
    
    col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
    col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
    col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
    
}

double cc(float orig)
{
    if(orig == 0)
    {
        return 255;
    }
    double off = orig*250/.5;
    if (orig>0)
    {
        return (off +250);
    }
    else
    {
        return (250 -off);
    }
    
    
}

double cy(float orig)
{
    if(orig == 0)
    {
        return 255;
    }
    double off = orig*250/.5;
    if (orig<0)
    {
        return (off*-1 +250);
    }
    else
    {
        return (250 -off);
    }
    
    std::cout <<  std::endl  << "helllo" << (250 -off) << std::endl << "test\n";
}

int main() {
    
    const int H = 500;
    const int W = 500;
    //convet to floats
    const Vec3 white(1.0, 1.0, 1.0);
    const Vec3 black(.1, .1, .1);
    const Vec3 red(1.0, 0, 0);
    const Vec3 green(0, .3, 0);
    const Vec3 ambient(.1, .1, .1);
    
    
    //const Sphere2 sphere2(Vec3(W*0.5 + 50, 450, 50), 50);
    //const Sphere sphere(Vec3(400,390, 50), 50);
    //const Triangle triangle(Vec3(100, 290, 60), Vec3(150, 0, 50), Vec3(250, 330, 160));
    
  
   
   // const Sphere light(Vec3(500, 250, 50), 1);
    
    const Sphere sphere(Vec3(.35, 0, -1) , .05);
    const Sphere2 sphere2(Vec3(.2), 0, -1), .75);
   const Triangle triangle(Vec3(.3,-.3,-.4), Vec3(0, .3, -.1), Vec3(cc(-.3), cy(-.3), cc(.2)));

    const Sphere light(Vec3(1, 0, 0), 1);
    
    std::ofstream out("out.ppm");
    out << "P3\n" << W << ' ' << H << ' ' << "255\n";
    
    double t;
    double t2;
    double tTriangle;
    Vec3 pix_col(black);
    
    Vec3 CameraLookAt = Vec3( 0, 0, 0);
    Vec3 CameraLookFrom = Vec3 (0, 0, 1);
    Vec3 diff = CameraLookAt - CameraLookFrom;
    
    
    double d = sqrt(diff.x*diff.x + diff.y * diff.y + diff.z * diff.z);
    
    double FieldOfView  = 28;
    
    //in radians !!!!!!
    double width = tan(FieldOfView)*d;
    
    double pixscale = 2*width/500;
    
    
    
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            pix_col = black;
            
            Vec3 pixelPoint = Vec3( (- width + (x *pixscale) + .5 * pixscale), ( width - (y *pixscale) - .5 * pixscale), 0);
            
            Vec3   rayDirection = (pixelPoint - CameraLookFrom).normalize();
            
            //0001 for first scene
            const Ray ray(Vec3(0,0,1),  rayDirection);
            
            
            
            
             if (sphere2.intersect(ray, t)) {
                const Vec3 pi = ray.o + ray.d*t;
                const Vec3 L = light.c - pi;
                const Vec3 N = sphere2.getNormal(pi);
                const double dt = dot(L.normalize(), N.normalize());
                //shadow ray check here
                
                /*
                 shadow ray direction (l – p) / ||l – p||
                 Distance to the light from the intersection point: tl = ||l – p||
                 
                 
                 In theory, the possible range of t-values is t  [0, tl]
                 Due to numerical (floating-point) error, test in the range t  [ε, t1] where ε  is some small value (such as 2–16)
                 
                 
                 */
                pix_col = (red + white*dt) * 0.5;
                clamp255(pix_col);
            }
            
            if (sphere.intersect(ray, t2)) {
                const Vec3 pi = ray.o + ray.d*t2;
                const Vec3 L = Vec3(1, 0, 0);
                const Vec3 N = sphere.getNormal(pi);
                const double dt = dot(L, N.normalize());
                
                //shadow ray check here~!!!!!!!!
                // ambient lighit
                pix_col = (red + white*dt ) * 0.5;
                clamp255(pix_col);
            }
            if (triangle.intersect(ray, tTriangle))
            {

                //assing two different t values
                const Vec3 pi = ray.o + ray.d*tTriangle;
                const Vec3 L = light.c - pi;
                const double dt = dot(L.normalize(), triangle.getNormal().normalize());

                //shadow raytrace
                pix_col = (green + white*dt) * 0.5;
                clamp255(pix_col);
            }
            
            out << (int)pix_col.x << ' '
            << (int)pix_col.y << ' '
            << (int)pix_col.z << '\n';
        }
    }
}


