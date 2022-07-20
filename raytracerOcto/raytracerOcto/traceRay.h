#pragma once


void traceSpheres(Scene& scene, vec3& s, vec3& v, vec3& ip, vec3& N, vec3& color, float& closestT, bool& inter)
{
    unsigned ci=(unsigned)-1;
    for(unsigned i=0;i<scene.spheres.size();++i){
        auto q = s - scene.spheres[i].c;
        auto A = dot(v,v);
        auto B = 2*dot(q,v);
        auto C = dot(q,q)-scene.spheres[i].r*scene.spheres[i].r;
        auto disc = B*B-4*A*C;
            if( disc < 0 )
            continue;
        disc = sqrt(disc);
        auto denom = 1.0 / (2.0*A);
        auto t1 = (-B + disc) * denom;
        auto t2 = (-B - disc) * denom;
        float t;
        if( t1 < 0 && t2 < 0 ){
            continue;
        } else if( t1 < 0 )
            t = t2;
        else if( t2 < 0 )
            t = t1;
        else
            t = ( (t1<t2) ? t1:t2 );
            
        if( t < closestT ){
            closestT = t;
            ci = i;
        }
    }
    
    if( ci == (unsigned)-1 ){
        inter=false;
        return;
    }
    ip = s + closestT * v;
    N = normalize(ip - scene.spheres[ci].c);
    color = scene.spheres[ci].color;
    inter=true;
}

void traceTriangles(Scene& scene, vec3& s, vec3& v, vec3& ipC, vec3& N, 
            vec3& color, float& closestT, bool& inter)
{
    for(auto& M : scene.meshes ){
        for(unsigned i=0;i<M.triangles.size();++i){
            Triangle& T = M.triangles[i];
            float denom = dot(T.N,v);
            if( denom == 0.0 )
                continue;
            float numer = -(T.D + dot(T.N,s) );
            float t = numer/denom;
            if( t < 0 )
                continue;
            auto ip = s + t * v;
            float A=0.0;
            for(unsigned j=0;j<3;++j){
                A += length( cross(T.e[j], ip-T.p[j] ) );
            }
            A *= T.oneOverTwiceArea;
            if( A > 1.001 )
                continue;
            if( t < closestT ){
                ipC = ip;
                N = T.N;
                color = M.color;
                closestT = t;
                inter = true;
            }
        }
    }
}

bool traceRay(Scene& scene, vec3& s, vec3& v, vec3& ip, vec3& N, vec3& color)
{
    float closestT = 1E99;
    bool inter;
    traceSpheres(scene,s,v,ip,N,color,closestT,inter);
    traceTriangles(scene,s,v,ip,N,color,closestT,inter);
    return inter;
}
