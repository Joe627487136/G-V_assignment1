#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
using namespace std;

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx( const Vector3f& lhs, const Vector3f& rhs )
    {
        const float eps = 1e-8f;
        return ( lhs - rhs ).absSquared() < eps;
    }

	const Matrix4f Bernstein_Matrix = Matrix4f(1,-3,3,-1,
                                               0,3,-6,3,
                                               0,0,3,-3,
                                               0,0,0,1);

	const Matrix4f Bernstein_Matrix_inverse = Matrix4f(1,1,1,1,
		                                               0,1.f/3,2.f/3,1,
													   0,0,1.f/3,1,
													   0,0,0,1);
    
    const Matrix4f B_Spline_Matrix = Matrix4f(1.f/6,-1.f/2,1.f/2,-1.f/6,
                                              2.f/3,0,-1,1.f/2,
                                              1.f/6,1.f/2,1.f/2,-1.f/2,
                                              0,0,0,1.f/6);

    
}
    

Curve evalBezier( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 || P.size() % 3 != 1 )
    {
		cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit( 0 );
    }
    
    Curve R;
    Vector3f B0 = Vector3f(0,0,1);

    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }
    
	unsigned num_of_steps;

    for (unsigned i=0; i < (P.size() - 1); i+=3){
        // Geometry matrix
		Vector3f P1 = P[i];
		Vector3f P2 = P[i+1];
		Vector3f P3 = P[i+2];
		Vector3f P4 = P[i+3];
		/*
		if (i>=P.size()-4){
			num_of_steps = steps+1;
		}else{
			num_of_steps = steps;
		}
		*/
        for (unsigned s=0; s<=steps; s++){
            float t = s*1/float(steps);
            Vector4f basis(1,t,t*t,t*t*t);
			Vector4f basis_d(0,1,2*t,3*t*t);
            // V = Q(t) = GBT(t)
			Vector4f V = Matrix4f(Vector4f(P1,0),
				                  Vector4f(P2,0),
                                  Vector4f(P3,0),
                                  Vector4f(P4,0))*Bernstein_Matrix*basis;
            // T = Q'(t).normalized()
			Vector4f T = Matrix4f(Vector4f(P1,0),
                                  Vector4f(P2,0),
                                  Vector4f(P3,0),
                                  Vector4f(P4,0))*Bernstein_Matrix*basis_d;
            Vector3f N;
            Vector3f B;


            
            if (i==0 && s==0){
				if (approx(B0, T.xyz())){
					B0 = Vector3f(1,0,0);
				}
			}

			CurvePoint curve_points;
			curve_points.V = V.xyz();
			curve_points.T = T.xyz().normalized();
			if (R.size()==0){
				curve_points.N = Vector3f::cross(B0,curve_points.T).normalized();
			}else{
				curve_points.N = Vector3f::cross(R[R.size()-1].B,curve_points.T).normalized();
			}
			curve_points.B = Vector3f::cross(curve_points.T,curve_points.N).normalized();
			R.push_back(curve_points);
		}
			
    }
	if (approx(R[0].V, R[R.size()-1].V) &&
		approx(R[0].T, R[R.size()-1].T) &&
		!approx(R[0].N, R[R.size()-1].N)){
			float angle = acos(Vector3f::dot(R[0].N, R[R.size()-1].N));
			for (unsigned i = 0; i < R.size(); i++){
				Matrix3f rotM = Matrix3f::rotation(R[i].T, angle * (.5-float(i)/R.size()));
				R[i].N = rotM*R[i].N;
				R[i].B = rotM*R[i].B;
			}
	}
    cerr << "\t>>> Steps (type steps): " << steps << endl;

    return R;
}

Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }

	Curve R;
	vector<Vector3f> Bernstein_Points;

	for( unsigned i = 0; i < P.size()-3; i++){
		// Geometry matrix
		Vector3f P1 = P[i];
		Vector3f P2 = P[i+1];
		Vector3f P3 = P[i+2];
		Vector3f P4 = P[i+3];
		Matrix4f bezierpts = Matrix4f(Vector4f(P1,0),
			                          Vector4f(P2,0),
			                          Vector4f(P3,0),
			                          Vector4f(P4,0))*B_Spline_Matrix*Bernstein_Matrix_inverse;
		Bernstein_Points.push_back(bezierpts.getCol(0).xyz());
		Bernstein_Points.push_back(bezierpts.getCol(1).xyz());
		Bernstein_Points.push_back(bezierpts.getCol(2).xyz());
		if (i == P.size()-4){
			Bernstein_Points.push_back(bezierpts.getCol(3).xyz());
		}
	}
	return evalBezier(Bernstein_Points, steps);
}

Curve evalCircle( float radius, unsigned steps )
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).
    
    // Preallocate a curve with steps+1 CurvePoints
    Curve R( steps+1 );

    // Fill it in counterclockwise
    for( unsigned i = 0; i <= steps; ++i )
    {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * float( i ) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f( cos(t), sin(t), 0 );
        
        // Tangent vector is first derivative
        R[i].T = Vector3f( -sin(t), cos(t), 0 );
        
        // Normal vector is second derivative
        R[i].N = Vector3f( -cos(t), -sin(t), 0 );

        // Finally, binormal is facing up.
        R[i].B = Vector3f( 0, 0, 1 );
    }

    return R;
}

void drawCurve( const Curve& curve, float framesize )
{
    // Save current state of OpenGL
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    // Setup for line drawing
    glDisable( GL_LIGHTING ); 
    glColor4f( 1, 1, 1, 1 );
    glLineWidth( 1 );
    
    // Draw curve
    glBegin( GL_LINE_STRIP );
    for( unsigned i = 0; i < curve.size(); ++i )
    {
        glVertex( curve[ i ].V );
    }
    glEnd();

    glLineWidth( 1 );

    // Draw coordinate frames if framesize nonzero
    if( framesize != 0.0f )
    {
        Matrix4f M;

        for( unsigned i = 0; i < curve.size(); ++i )
        {
            M.setCol( 0, Vector4f( curve[i].N, 0 ) );
            M.setCol( 1, Vector4f( curve[i].B, 0 ) );
            M.setCol( 2, Vector4f( curve[i].T, 0 ) );
            M.setCol( 3, Vector4f( curve[i].V, 1 ) );

            glPushMatrix();
            glMultMatrixf( M );
            glScaled( framesize, framesize, framesize );
            glBegin( GL_LINES );
            glColor3f( 1, 0, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 1, 0, 0 );
            glColor3f( 0, 1, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 1, 0 );
            glColor3f( 0, 0, 1 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 0, 1 );
            glEnd();
            glPopMatrix();
        }
    }
    
    // Pop state
    glPopAttrib();
}

