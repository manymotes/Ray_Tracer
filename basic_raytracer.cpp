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
     v1.x * v2.y - v1.z * v2.x
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


//struct Sphere {
//    Vec3 c;
//    double r;
//    Sphere(const Vec3& c, double r) : c(c), r(r) {}
//    Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
//
//
//    bool sIntersect(const Ray& ray, double &t) const {
//        const Vec3 o = ray.o;
//        const Vec3 d = ray.d;
//        const Vec3 oc = o - c;
//        const double b = 2 * dot(oc, d);
//        const double c = dot(oc, oc) - r*r;
//        double disc = b*b - 4 * c;
//        if (disc < 1e-4) return false;
//        disc = sqrt(disc);
//        const double t0 = -b - disc;
//        const double t1 = -b + disc;
//
//        const double tca = dot(d.normalize(), oc);
//        const double thc = r*r - oc.findLength() + tca * tca;
//        // t = tca - thc;
//        //t = (t0 < t1) ? t0 : t1;
//        return true;
//    }
//
//
//
//    bool intersect(const Ray& ray, double &t) const {
//        const Vec3 o = ray.o;
//        const Vec3 d = ray.d;
//        const Vec3 oc = o - c;
//        const double b = 2 * dot(oc, d);
//        const double c = dot(oc, oc) - r*r;
//        double disc = b*b - 4 * c;
//        if (disc < 1e-4) return false;
//        disc = sqrt(disc);
//        const double t0 = -b - disc;
//        const double t1 = -b + disc;
//
//        const double tca = dot(d.normalize(), oc);
//        const double thc = r*r - oc.findLength() + tca * tca;
//       // t = tca - thc;
//        t = (t0 < t1) ? t0 : t1;
//        return true;
//    }
//};


struct Sphere2 {
    Vec3 c;
    double r;
    Sphere2(const Vec3& c, double r) : c(c), r(r) {}
    Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
   
    bool sIntersect(const Ray& ray, double &t) const {
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
        //t = (t0 < t1) ? t0 : t1;
        return true;
    }
    
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
    
    bool sIntersect(const Ray& ray, double &t_out) const{
        const Vec3 o = ray.o;
        const Vec3 d = ray.d;
       
        return false;
        
    }
    
    
    
    bool intersect(const Ray& ray, double &t_out) const {
        
//
//        Vec3 pn = getNormal();
//
//        const Vec3 o = ray.o;
//        const Vec3 d = ray.d;
//        //Vec3 normal = CrossProduct(b-a,tA),(b-c)).normalize();
//        double dis =  -dot(pn, b);
//
//        double t = -((pn.x * o.x) + (pn.y * o.y) + (pn.z * o.z) + dis)/( (pn.x * d.x) + (pn.y * d.y) + (pn.z * d.z));
//
//        Vec3 intersect = (o + MultiplyScalar(d, t));
//
//        bool signA = sign(dot(CrossProduct((b - a),(intersect - a)),pn));
//        bool signB = sign(dot(CrossProduct((c - b),(intersect - b)),pn));
//        bool signC = sign(dot(CrossProduct((a - c),(intersect - c)),pn));
//
//
//        //help with this t value
//        //t_out = (t0 < t1) ? t0 : t1;
//        t_out = t;
//        return ((signA == signB) && (signA == signC) && t>0.1);
//
//
//
        
        

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
//        if (ri != Vec2(0,0))
//        {
//            throw noexcept;
//        }

        //u and v prime are g0 g2 and g2
        //singholde, look at the v coordinate of the first point. which is 3 in the example.
        // sign holder = 1, becuase 3 is positive.
        
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
    //convet to floats
    const Vec3 white(1, 1, 1);
    const Vec3 black(.1, .1, .1);
    const Vec3 red(1, 0, 0);
    const Vec3 green(0, .5, 0);
    const Vec3 ambient(.1, .1, .1);
    const Vec3 blue(0, 0, 1);
    
//    const Vec3 white(255, 255, 255);
//    const Vec3 black(.1, .1, .1);
//    const Vec3 red(255, 0, 0);
//    const Vec3 green(0, 300, 0);
//    const Vec3 ambient(.1, .1, .1);
    
    
    
    //const Sphere2 sphere2(Vec3(W*0.5 + 50, 450, 50), 50);
    //const Sphere sphere(Vec3(400,390, 50), 50);
    //const Triangle triangle(Vec3(100, 290, 60), Vec3(150, 0, 50), Vec3(250, 330, 160));
    
  
   
   // const Sphere light(Vec3(500, 250, 50), 1);
    
    const Sphere sphere(Vec3(.35, 0, -.1) , .05);
    const Sphere2 sphere2(Vec3(.2, 0, -.1), .075);
   const Triangle triangle(Vec3(.3,-.3,-.4), Vec3(0, .3, -.1), Vec3(-.3, -.3, .2));

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
    
    double FieldOfView = 28;
    
    //in radians !!!!!!
    double width = tan(FieldOfView * 3.14159265 / 180)*d;
    
    double pixscale = 2*width/500;
    
    
    
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            pix_col = black;
            
            Vec3 pixelPoint = Vec3( (- width + (x *pixscale) + .5 * pixscale), ( width - (y *pixscale) - .5 * pixscale), 0);
            
            Vec3   rayDirection = (pixelPoint - CameraLookFrom).normalize();
            
            //0001 for first scene
            const Ray ray(Vec3(0,0,1),  rayDirection);
            
            
            
            
             if (sphere2.intersect(ray, t2)) {
                const Vec3 pi = ray.o + ray.d*t2;
                const Vec3 L = Vec3(1, 0, 0);
                const Vec3 N = sphere2.getNormal(pi).normalize();
                const double dt = dot(L.normalize(), N.normalize());
                 
                 
                //shadow ray check here
                
                 
                /*
                 shadow ray direction (l – p) / ||l – p||
                 Distance to the light from the intersection point: tl = ||l – p||
                 
                 
                 In theory, the possible range of t-values is t  [0, tl]
                 Due to numerical (floating-point) error, test in the range t  [ε, t1] where ε  is some small value (such as 2–16)
                 
                 
                 */
                 
                 Ray shadowRay(pi, (    (L - pi)   /   ((L - pi).normalize())  ));
                 
            // if (!sphere.sIntersect(shadowRay, t2))
          //   {
                     
                pix_col = (red + white*dt) * 0.5;
                scale255(pix_col);
              //   }
             
            }
            
            if (sphere.intersect(ray, t)) {
                
                
                //shadow ray check here~!!!!!!!!
                // ambient lighit
                
                /*
                 cr object color
                 ca ambient
                 cl light source
                 n normal
                 l direction to light from point
                 cp specular
                 e direction to eye
                 r reflection from light bout normal
                 p phong constant
                 */
                
                
                const Vec3 pi = ray.o + ray.d*t;
                const Vec3 L = Vec3(1, 0, 0);
                const Vec3 N = sphere.getNormal(pi).normalize();
                const double dt = dot(L, N.normalize());
                const Vec3 diffuse = Vec3(1, 1, 1);
                const Vec3 specHigh = Vec3(1,1,1);
                
                
                Vec3 r =  (MultiplyScalar(N,2) * dot(N, L)  - L).normalize();
                //float spec = pow( fmax( 0.0, dot( CameraLookFrom, r)), 32 );
               
                
              //legit working on
                Vec3 v = (CameraLookFrom - pi).normalize();
                
                pix_col = white * (ambient + diffuse * fmax(0.0, dot(N,L)) + specHigh * white * fmax(0.0, dot(v,r)))  ;
                //pix_col = (ambient * red + diffuse * red * dot(N, L) + specHigh * white * fmax(0.0, dot(v,r)));
                scale255(pix_col);
                
                
               // pix_col = CrossProduct(red, ()
                
               // pix_col =(red + white*dt ) * 0.5;
                //clamp255(pix_col);
                
            }
            
            if (triangle.intersect(ray, tTriangle) )
            {

                //assing two different t values
                const Vec3 pi = ray.o + ray.d*tTriangle;
                const Vec3 L = Vec3(1, 0, 0);
                
                Ray shadowRay(pi, (    (L - pi)   /   ((L - pi).normalize())  ));
                
                
                
                //???????????????
               //if (!sphere2.sIntersect(shadowRay, t2)){
                const double dt = dot(L.normalize(), triangle.getNormal().normalize());
                const Vec3 N = triangle.getNormal().normalize();
                
                
                Vec3 r = MultiplyScalar( MultiplyScalar(N, 2), dot(N, L) ) - L;
                //float spec = pow( fmax( 0.0, dot( CameraLookFrom, r)), 32 );
                Vec3 specHigh = Vec3(1,1,1);
                
                
                const Vec3 diffuse = Vec3(1, 1, 1);
                Vec3 v = (CameraLookFrom - pi).normalize();
                
             //  pix_col = (ambient * diffuse + diffuse * white * abs(dot(N, L)) * blue + specHigh * white * pow( abs(dot(CameraLookFrom,    (( L * 2 *  abs(dot( N, L))   )  - L ).normalize()  )   ),  32)   );
                 
                pix_col =  white* (ambient + diffuse * fmax( 0.0, dot(N,L)));
                
                scale255(pix_col);
                
                
                
                //shadow raytrace
            //    pix_col = white;//(green + white*dt) * 0.5;
         //     scale255(pix_col);
               
            //shadow if statement end
            //}
            }
            
            out << (int)pix_col.x << ' '
            << (int)pix_col.y << ' '
            << (int)pix_col.z << '\n';
        }
    }
}


