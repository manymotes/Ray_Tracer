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
     Vec3 operator * (Vec3 v) const { return Vec3(x*v.x, y*v.y, z*v.z); }
    Vec3 operator / (double d) const { return Vec3(x/d, y/d, z/d); }
    Vec3 operator / (Vec3 v) const { return Vec3(x/v.x, y/v.y, z/v.z); }
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

struct Vec2
{
    Vec2(double x, double y) : x(x), y(y) {}
    float x;
    float y;
    
    Vec2 operator - (const Vec2& v) const { return Vec2(x-v.x, y-v.y); }
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
     v1.x * v2.y - v1.y * v2.x
     );
    
    return toReturn;
}

struct Ray {
    Vec3 o,d;
    Ray(const Vec3& o, const Vec3& d) : o(o), d(d) {}
};

struct Sphere {
    Vec3 c;
    double r;
    Sphere(const Vec3& c, double r) : c(c), r(r) {}
    Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
    bool intersect(const Ray& ray, double &t) const {
        const Vec3 o = ray.o;
        const Vec3 d = ray.d;
        const Vec3 oc = c - o;
        const double tca = dot(d.normalize(), oc);
        if (tca < 0)
        {
            return false;
        }
        const double thc = r*r - oc.findLength()*oc.findLength() + tca * tca;
        if (thc < 0)
        {
            return false;
        }
        t = tca - sqrt(thc);

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
        const Vec3 oc = c - o;
        const double tca = dot(d.normalize(), oc);
        if (tca < 0)
        {
            return false;
        }
        const double thc = r*r - oc.findLength()*oc.findLength() + tca * tca;
        if (thc < 0)
        {
            return false;
        }
        t = tca - sqrt(thc);
        return true;
    }
};

struct Sphere3 {
    Vec3 c;
    double r;
    Sphere3(const Vec3& c, double r) : c(c), r(r) {}
    Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
    
    bool intersect(const Ray& ray, double &t) const {
        const Vec3 o = ray.o;
        const Vec3 d = ray.d;
        const Vec3 oc = c - o;
        const double tca = dot(d.normalize(), oc);
        if (tca < 0)
        {
            return false;
        }
        const double thc = r*r - oc.findLength()*oc.findLength() + tca * tca;
        if (thc < 0)
        {
            return false;
        }
        t = tca - sqrt(thc);
        
        return true;
    }
    
    
};

struct Triangle {
    Vec3 a;
    Vec3 b;
    Vec3 c;
    
 
    Triangle(const Vec3& a, const Vec3& b, const Vec3& c) : a(a), b(b), c(c)  {
    }

    Vec3 getNormal() const
    {
        Vec3 v1 = b-a;
        Vec3 v2 = c-a;
        Vec3 pn = CrossProduct(v1, v2);
        return pn.normalize();
    }
 
    bool intersect(const Ray& ray, double &t_out) const {
        
 
        const Vec3 o = ray.o;
        const Vec3 d = ray.d;
        //        const Vec3 oc = o - c;
        //        const double b = 2 * dot(oc, d);
        //        const double c = dot(oc, oc) - r*r;
        //        double disc = b*b - 4 * c;
//        Vec3 v1 = a-b;
//        Vec3 v2 = b-c;
        Vec3 pn = getNormal();
        
        ////!!!!! find dominant axis!!!!!
        std::string dominatAxis;
        if (abs(pn.x) > abs(pn.y))
        {
            if(abs(pn.x) >= abs(pn.z))
            {
                dominatAxis = "x";
            }
            else{
                dominatAxis ="z";
            }
        }
        else{
            if (abs(pn.y) >= abs(pn.z))
            {
                dominatAxis = "y";
            }
            else
            {
                dominatAxis = "z";
            }
        }

        //gets intersection of ray and plane

        double dis = -dot(pn, b);
        double t = -((pn.x * o.x) + (pn.y * o.y) + (pn.z * o.z) + dis)/( (pn.x * d.x) + (pn.y * d.y) + (pn.z * d.z));
        t_out = t;
        
        Vec2 projectedVertices0 = Vec2(a.y, a.z);
        Vec2 projectedVertices1 =Vec2(b.y, b.z);
        Vec2 projectedVertices2 = Vec2(c.y, c.z);
        Vec3 intersect = Vec3((o.x + d.x*t), (o.y + d.y*t), (o.z + d.z*t));
        Vec2 ri = Vec2(intersect.y, intersect.z);
        
        if (dominatAxis == "x")
        {
            projectedVertices0 = Vec2(a.y, a.z);
            projectedVertices1 =Vec2(b.y, b.z);
            projectedVertices2 = Vec2(c.y, c.z);
            ri = Vec2(intersect.y, intersect.z);
        }
        else if (dominatAxis == "y")
        {
            projectedVertices0 = Vec2(a.x, a.z);
            projectedVertices1 = Vec2(b.x, b.z);
            projectedVertices2 = Vec2(c.x, c.z);
            ri = Vec2(intersect.x, intersect.z);
        }
        else if(dominatAxis == "z")
        {
            projectedVertices0 = Vec2(a.x, a.y);
            projectedVertices1 = Vec2(b.x, b.y);
            projectedVertices2 = Vec2(c.x, c.y);
            ri = Vec2(intersect.x, intersect.y);
        }


        //translate verticies by -ri
        projectedVertices0 = projectedVertices0 - ri;
        projectedVertices1 = projectedVertices1 - ri;
        projectedVertices2 = projectedVertices2 - ri;
        
        Vec2 trianlgeEdges[4] = {projectedVertices0, projectedVertices1, projectedVertices2, projectedVertices0};
        
        ri = ri - ri;
        
        int signHolder;
        if (trianlgeEdges[0].y < 0 )
        {
            signHolder = -1;
        }
        else
        {
            signHolder = 1;
        }
        
        int nextsing;
        int numCrossing = 0;
        
        for( int j = 0; j < 3; j++)
        {
            if (trianlgeEdges[j+1].y < 0)
            {
                nextsing = -1;
            }
            else {
                nextsing = 1;
            }
            
            if (signHolder != nextsing)
            {
                if(trianlgeEdges[j].x > 0 && trianlgeEdges[j+1].x > 0)
                {
                    numCrossing +=1;
                }
                else if(trianlgeEdges[j].x > 0 || trianlgeEdges[j+1].x > 0)
                {
                    float uCross = trianlgeEdges[j].x - trianlgeEdges[j].y * (trianlgeEdges[j+1].x - trianlgeEdges[j].x)/(trianlgeEdges[j+1].y - trianlgeEdges[j].y);
                    
                    if (uCross > 0)
                    {
                        numCrossing += 1;
                    }
                }
            }
            signHolder = nextsing;
        }
        
        if (numCrossing % 2== 0)
        {
            return false;
        }
        return true;
    }
    

};

struct Triangle2 {
    Vec3 a;
    Vec3 b;
    Vec3 c;
    
    
   Vec3 getLocationWithMagnitude(Ray ray, double magnitude) const
    {
        return ray.o  + (ray.d * magnitude);
    }
    
//     Vec3 getPlaneIntersectionVector(Ray ray) const
//    {
//        Vec3 zero = Vec3(0,0,0);
//        Vec3 N = getNormal();
//        double epsilon = 0.00000001;
//        Vec3 w0 = ray.o - b;
//        double numerator = -dot(N, w0);
//        double denominator = dot(N, ray.d);
//
//        if (abs(denominator) < epsilon)
//        {
//            //ray lies in triangle plane
//            if (numerator == 0)
//            {
//                return zero;
//            }
//            //ray is disjoint from plane
//            else
//            {
//                return zero;
//            }
//        }
//
//        double intersectionDistance = numerator / denominator;
//         return (intersectionDistance >= 0) ? getLocationWithMagnitude(ray, intersectionDistance) : zero;
//
//    }
    
    Triangle2(const Vec3& a, const Vec3& b, const Vec3& c) : a(a), b(b), c(c)  {
    }
    
    Vec3 getNormal() const
    {
        Vec3 v1 = b-a;
        Vec3 v2 = c-a;
        Vec3 pn = CrossProduct(v1, v2);
        return pn.normalize();
    }

    bool intersect(const Ray& ray, double &t_out) const {
        
    
        
        const Vec3 o = ray.o;
        const Vec3 d = ray.d;
        //        const Vec3 oc = o - c;
        //        const double b = 2 * dot(oc, d);
        //        const double c = dot(oc, oc) - r*r;
        //        double disc = b*b - 4 * c;
        //        Vec3 v1 = a-b;
        //        Vec3 v2 = b-c;
        Vec3 pn = getNormal();
        
        ////!!!!! find dominant axis!!!!!
        std::string dominatAxis;
        if (abs(pn.x) > abs(pn.y))
        {
            if(abs(pn.x) >= abs(pn.z))
            {
                dominatAxis = "x";
            }
            else{
                dominatAxis ="z";
            }
        }
        else{
            if (abs(pn.y) >= abs(pn.z))
            {
                dominatAxis = "y";
            }
            else
            {
                dominatAxis = "z";
            }
        }
        
        //gets intersection of ray and plane
        
        double dis = -dot(pn, b);
        double t = -((pn.x * o.x) + (pn.y * o.y) + (pn.z * o.z) + dis)/( (pn.x * d.x) + (pn.y * d.y) + (pn.z * d.z));
        t_out = t;
        
        Vec2 projectedVertices0 = Vec2(a.y, a.z);
        Vec2 projectedVertices1 =Vec2(b.y, b.z);
        Vec2 projectedVertices2 = Vec2(c.y, c.z);
        Vec3 intersect = Vec3((o.x + d.x*t), (o.y + d.y*t), (o.z + d.z*t));
        Vec2 ri = Vec2(intersect.y, intersect.z);
        
        if (dominatAxis == "x")
        {
            projectedVertices0 = Vec2(a.y, a.z);
            projectedVertices1 =Vec2(b.y, b.z);
            projectedVertices2 = Vec2(c.y, c.z);
            ri = Vec2(intersect.y, intersect.z);
        }
        else if (dominatAxis == "y")
        {
            projectedVertices0 = Vec2(a.x, a.z);
            projectedVertices1 = Vec2(b.x, b.z);
            projectedVertices2 = Vec2(c.x, c.z);
            ri = Vec2(intersect.x, intersect.z);
        }
        else if(dominatAxis == "z")
        {
            projectedVertices0 = Vec2(a.x, a.y);
            projectedVertices1 = Vec2(b.x, b.y);
            projectedVertices2 = Vec2(c.x, c.y);
            ri = Vec2(intersect.x, intersect.y);
        }
        
        
        //translate verticies by -ri
        projectedVertices0 = projectedVertices0 - ri;
        projectedVertices1 = projectedVertices1 - ri;
        projectedVertices2 = projectedVertices2 - ri;
        
        Vec2 trianlgeEdges[4] = {projectedVertices0, projectedVertices1, projectedVertices2, projectedVertices0};
        
        ri = ri - ri;
        
        int signHolder;
        if (trianlgeEdges[0].y < 0 )
        {
            signHolder = -1;
        }
        else
        {
            signHolder = 1;
        }
        
        int nextsing;
        int numCrossing = 0;
        
        for( int j = 0; j < 3; j++)
        {
            if (trianlgeEdges[j+1].y < 0)
            {
                nextsing = -1;
            }
            else {
                nextsing = 1;
            }
            
            if (signHolder != nextsing)
            {
                if(trianlgeEdges[j].x > 0 && trianlgeEdges[j+1].x > 0)
                {
                    numCrossing +=1;
                }
                else if(trianlgeEdges[j].x > 0 || trianlgeEdges[j+1].x > 0)
                {
                    float uCross = trianlgeEdges[j].x - trianlgeEdges[j].y * (trianlgeEdges[j+1].x - trianlgeEdges[j].x)/(trianlgeEdges[j+1].y - trianlgeEdges[j].y);
                    
                    if (uCross > 0)
                    {
                        numCrossing += 1;
                    }
                }
            }
            signHolder = nextsing;
        }
        
        if (numCrossing % 2== 0)
        {
            return false;
        }
        return true;
    }
    
    
};


void clamp255(Vec3& col) {
    col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
    col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
    col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
    
}

void scale255(Vec3& col)
{
    col.x = col.x * 255;
    col.y = col.y * 255;
    col.z = col.z * 255;
    
    clamp255(col);
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
   
    const Vec3 white(1, 1, 1);
    const Vec3 black(.2, .2, .2);
    const Vec3 red(1, 0, 0);
    const Vec3 green(0, 1, 0);
    const Vec3 ambient(.1, .1, .1);
    const Vec3 blue(0, 0, 1);
    const Vec3 yellow(1, 1, 0);
    const Vec3 dark(.1, .1, .1);
    
    

    const Sphere sphere(Vec3(.35, 0, -.1) , .05);
    const Sphere2 sphere2(Vec3(.2, 0, -.1), .075);
    const Sphere3 sphere3(Vec3(-.6, 0, 0), .3);
    const Triangle triangle(Vec3(.3,-.3,-.4), Vec3(0, .3, -.1), Vec3(-.3, -.3, .2));
    const Triangle2 triangle2(Vec3(-.2,.1,.1), Vec3(-.2, -.5, .2), Vec3(-.2, .1, -.3));
    
    
    const Sphere light(Vec3(1, 0, 0), 1);
    
    std::ofstream out("out.ppm");
    out << "P3\n" << W << ' ' << H << ' ' << "255\n";
    
    double t;
    double t2;
    double t3;
    double tTriangle;
    double tTriangle2;
    Vec3 pix_col(black);
    
    Vec3 CameraLookAt = Vec3( 0, 0, 0);
    Vec3 CameraLookFrom = Vec3 (0, 0, 1);
    Vec3 diff = CameraLookAt - CameraLookFrom;
    
    
    double d = sqrt(diff.x*diff.x + diff.y * diff.y + diff.z * diff.z);
    
    double FieldOfView = 28;
    
    //in radians !!!!!!
    double width = tan(FieldOfView * 3.14159265 / 180)*d;
    
    double pixscale = 2*width/500;
    
    
    
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            pix_col = black;
            scale255(pix_col);
            
            Vec3 pixelPoint = Vec3( (- width + (x *pixscale) + .5 * pixscale), ( width - (y *pixscale) - .5 * pixscale), 0);
            
            Vec3   rayDirection = (pixelPoint - CameraLookFrom).normalize();
            
            //0001 for first scene
            const Ray ray(Vec3(0,0,1),  rayDirection);
            
            
            if (sphere.intersect(ray, t)) {
                
                
                const Vec3 pi = ray.o + ray.d*t;
                const Vec3 L = Vec3(1, 0, 0);
                const Vec3 N = sphere.getNormal(pi).normalize();
                const double dt = dot(L, N.normalize());
                const Vec3 diffuse = Vec3(1, 1, 1);
                const Vec3 specHigh = Vec3(1,1,1);
                const double phong = 4.0;
                
                
                Vec3 r =  (MultiplyScalar(N,2) * dot(N, L)  - L).normalize();
                
                Vec3 v = (CameraLookFrom - pi).normalize();
                
                pix_col = white * (ambient + diffuse * fmax(0.0, dot(N,L))) + specHigh * white * pow(fmax(0.0, dot(v,r)), phong);
                
                scale255(pix_col);
            }
            
            if (triangle.intersect(ray, tTriangle) )
            {
                
                //assing two different t values
                const Vec3 pi = ray.o + ray.d*tTriangle;
                const Vec3 L = Vec3(1, 0, 0);
                const double dt = dot(L.normalize(), triangle.getNormal().normalize());
                const Vec3 N = triangle.getNormal().normalize();
                
                
                Vec3 r = MultiplyScalar( MultiplyScalar(N, 2), dot(N, L) ) - L;
                
                Vec3 specHigh = Vec3(1,1,1);
                const double phong = 32;
                const Vec3 diffuse = Vec3(1, 1, 1);
                
                Vec3 v = (CameraLookFrom - pi).normalize();
                
                double bias = .0001;
                
                Ray shadowRay((pi + N * bias), (L - pi).normalize());
                if (sphere2.intersect(shadowRay, t2))
                {
                    pix_col = ambient * blue;
                }
                else
                {
                    pix_col = blue * (ambient + diffuse * fmax(0.0, dot(N,L))) + specHigh * white * pow(fmax(0.0, dot(v,r)), phong);
                }
                
                scale255(pix_col);
                
            }
            
             if (sphere2.intersect(ray, t2)) {
                const Vec3 pi = ray.o + ray.d*t2;
                const Vec3 L = Vec3(1, 0, 0);
                const Vec3 N = sphere2.getNormal(pi).normalize();
                const double dt = dot(L.normalize(), N.normalize());
                const Vec3 diffuse  = Vec3(1, 0, 0);
                 const Vec3 specHigh = Vec3(.5, 1, .5);
                 const double phong = 32;
                 
                 double bias = .0001;
                
                Ray shadowRay(pi + N * bias, (L - pi).normalize());
             if (sphere.intersect(shadowRay, t))
             {
                 pix_col = ambient * red;
             }
             else
             {
                 Vec3 v = (CameraLookFrom - pi).normalize();
                 Vec3 r =  (MultiplyScalar(N,2) * dot(N, L)  - L).normalize();
                 
                 pix_col = red * (ambient + diffuse * fmax(0.0, dot(N,L))) +  specHigh * white * pow(fmax(0.0, dot(v,r)), phong);
             }
                 
                 scale255(pix_col);
            }
            
          
            
            
            
            
            
            if (sphere3.intersect(ray, t3)) {
                const Vec3 pi = ray.o + ray.d*t3;
                const Vec3 L = Vec3(1, 0, 0);
                 Vec3 N = sphere3.getNormal(pi).normalize();
                const double dt = dot(L.normalize(), N.normalize());
                
               
                const Vec3 diffuse  = Vec3(0, 1, 0);
                const Vec3 specHigh = Vec3(.5, 1, .5);
                const double phong = 32;
                
                
                
                
                double bias = .0001;
                
                Ray shadowRay(pi + N * bias, (L - pi).normalize());
                
                
                if (triangle.intersect(shadowRay, tTriangle) || triangle2.intersect(shadowRay, tTriangle2))
                {
                    pix_col = ambient * green;
                }
                else
                {
               
                    Vec3 v = (CameraLookFrom - pi).normalize();
                    Vec3 r =  (MultiplyScalar(N,2) * dot(N, L)  - L).normalize();
                    
                    pix_col = green * (ambient + diffuse * fmax(0.0, dot(N,L))) + specHigh * white * pow(fmax(0.0, dot(v,r)), phong);
                    
                }
                    scale255(pix_col);
                
                
            }

            


            
            
            if (triangle2.intersect(ray, tTriangle2) )
            {
                
                //assing two different t values
                const Vec3 pi = ray.o + ray.d*tTriangle2;
                const Vec3 L = Vec3(1, 0, 0);
                
             
                //if (!sphere2.sIntersect(shadowRay, t2)){
                const double dt = dot(L.normalize(), triangle2.getNormal().normalize());
                const Vec3 N = triangle2.getNormal().normalize();
                
                
                Vec3 r = MultiplyScalar( MultiplyScalar(N, 2), dot(N, L) ) - L;
                
                Vec3 specHigh = Vec3(1,1,1);
                const Vec3 diffuse = Vec3(1, 1, 0);
                Vec3 v = (CameraLookFrom - pi).normalize();
                const double phong = 4;
                
                double bias = .0001;
                
                double distance = L.findLength();
                Ray shadowRay(pi + N * .002, ((L - pi).normalize())/distance);
                if (sphere2.intersect(shadowRay, t2) || triangle.intersect(shadowRay, tTriangle) )
                {
                    pix_col = ambient * yellow;
                }
                else
                {
                pix_col = yellow * (ambient + diffuse * fmax(0.0, dot(N,L))) + specHigh * white * pow(fmax(0.0, dot(v,r)), phong);
                }
                
                scale255(pix_col);
                
            }
            out << (int)pix_col.x << ' '
            << (int)pix_col.y << ' '
            << (int)pix_col.z << '\n';
        }
    }
}

/*
 t0 < 0
 */
