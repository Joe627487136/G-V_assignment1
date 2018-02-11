#include "surf.h"
#include "extra.h"
using namespace std;

namespace
{
    
    // We're only implenting swept surfaces where the profile curve is
    // flat on the xy-plane.  This is a check function.
    static bool checkFlat(const Curve &profile)
    {
        for (unsigned i=0; i<profile.size(); i++)
            if (profile[i].V[2] != 0.0 ||
                profile[i].T[2] != 0.0 ||
                profile[i].N[2] != 0.0)
                return false;
    
        return true;
    }
}

Surface makeSurfRev(const Curve &profile, unsigned steps)
{
    Surface surface;
    
    if (!checkFlat(profile))
    {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }
	for (unsigned step = 0; step < steps; step++) {
		float rotation = 2.0*M_PI*step/steps;
		Matrix4f M = Matrix4f::rotateY(rotation);
		Matrix3f transform = M.getSubmatrix3x3(0,0).inverse().transposed();

		for (unsigned i = 0; i < profile.size(); i++){

			CurvePoint curve_points = profile[i];
			Vector3f vertex = (M*Vector4f(curve_points.V,1)).xyz();
			Vector3f norm_vec = (-1*transform*curve_points.N).xyz().normalized();
			surface.VV.push_back(vertex);
			surface.VN.push_back(norm_vec);

			if (i < profile.size() - 1) {
				unsigned triangle_1 = i+step*profile.size();
				unsigned triangle_2 = (triangle_1+profile.size())%(steps*profile.size());
				Tup3u tup1(triangle_1, triangle_1+1, triangle_2);
				surface.VF.push_back(tup1);
				Tup3u tup2(triangle_1+1, triangle_2+1, triangle_2);
				surface.VF.push_back(tup2);
			}

		}
	}

    cerr << "\t>>> makeSurfRev called (but not implemented).\n\t>>> Returning empty surface." << endl;
 
    return surface;
}

Surface makeGenCyl(const Curve &profile, const Curve &sweep )
{
    Surface surface;

    if (!checkFlat(profile))
    {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // TODO: Here you should build the surface.  See surf.h for details.
	for (unsigned step = 0; step < sweep.size(); step++) {
		Matrix4f basis = Matrix4f(Vector4f(sweep[step].N, 0.0)
			                    , Vector4f(sweep[step].B, 0.0)
			                    , Vector4f(sweep[step].T, 0.0)
			                    , Vector4f(sweep[step].V, 1.0), true);
		Matrix3f norm = Matrix3f(basis.getCol(0).xyz(),basis.getCol(1).xyz(),basis.getCol(2).xyz(),true).inverse().transposed();
		for (unsigned i = 0; i < profile.size(); i++){
			Vector3f vertex = (basis*Vector4f(profile[i].V, 1.f)).xyz();
			Vector3f norm_vec = -1.0 * (norm * profile[i].N);
			surface.VV.push_back(vertex);
			surface.VN.push_back(norm_vec);
			if (i < profile.size() - 1) {
				unsigned triangle_1 = i + step * profile.size();
				unsigned triangle_2 = (triangle_1+profile.size())%(sweep.size()*profile.size());
				Tup3u tup1(triangle_1, triangle_1+1, triangle_2);
				surface.VF.push_back(tup1);
				Tup3u tup2(triangle_1+1, triangle_2+1, triangle_2);
				surface.VF.push_back(tup2);
			}
		}
	}

    return surface;
}

void drawSurface(const Surface &surface, bool shaded)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    if (shaded)
    {
        // This will use the current material color and light
        // positions.  Just set these in drawScene();
        glEnable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // This tells openGL to *not* draw backwards-facing triangles.
        // This is more efficient, and in addition it will help you
        // make sure that your triangles are drawn in the right order.
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
    }
    else
    {        
        glDisable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
        glColor4f(0.4f,0.4f,0.4f,1.f);
        glLineWidth(1);
    }

    glBegin(GL_TRIANGLES);
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        glNormal(surface.VN[surface.VF[i][0]]);
        glVertex(surface.VV[surface.VF[i][0]]);
        glNormal(surface.VN[surface.VF[i][1]]);
        glVertex(surface.VV[surface.VF[i][1]]);
        glNormal(surface.VN[surface.VF[i][2]]);
        glVertex(surface.VV[surface.VF[i][2]]);
    }
    glEnd();

    glPopAttrib();
}

void drawNormals(const Surface &surface, float len)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDisable(GL_LIGHTING);
    glColor4f(0,1,1,1);
    glLineWidth(1);

    glBegin(GL_LINES);
    for (unsigned i=0; i<surface.VV.size(); i++)
    {
        glVertex(surface.VV[i]);
        glVertex(surface.VV[i] + surface.VN[i] * len);
    }
    glEnd();

    glPopAttrib();
}

void outputObjFile(ostream &out, const Surface &surface)
{
    
    for (unsigned i=0; i<surface.VV.size(); i++)
        out << "v  "
            << surface.VV[i][0] << " "
            << surface.VV[i][1] << " "
            << surface.VV[i][2] << endl;

    for (unsigned i=0; i<surface.VN.size(); i++)
        out << "vn "
            << surface.VN[i][0] << " "
            << surface.VN[i][1] << " "
            << surface.VN[i][2] << endl;

    out << "vt  0 0 0" << endl;
    
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        out << "f  ";
        for (unsigned j=0; j<3; j++)
        {
            unsigned a = surface.VF[i][j]+1;
            out << a << "/" << "1" << "/" << a << " ";
        }
        out << endl;
    }
}
