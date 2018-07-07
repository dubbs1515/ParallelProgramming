typedef float4 point;
typedef float4 vector;
typedef float4 color;
typedef float4 sphere;


vector
Bounce( vector in, vector n )
{
	vector out = in - n*(vector)( 2.*dot(in.xyz, n.xyz) );
	out.w = 0.;
	return out;
}

vector
BounceSphere( point p, vector v, sphere s )
{
	vector n;
	n.xyz = fast_normalize( p.xyz - s.xyz );
	n.w = 0.;
	return Bounce( v, n );
}

bool
IsInsideSphere( point p, sphere s )
{
	float r = fast_length( p.xyz - s.xyz );
	return  ( r < s.w );
}

kernel
void
Particle( global point *dPobj, global vector *dVel, global color *dCobj )
{
	const float4 G       = (float4) ( 0., -40.8, -9.0, 0.0 ); // Manipulated gravity for effect
	const float  DT      = 0.08;	// Reduced time step
	
	const sphere Sphere1 = (sphere)( -0., -0., 0.,  1000. ); // Outer Sphere
	const sphere Sphere2 = (sphere)(-0., -0., 0., 300.);	// Inner Sphere

	int gid = get_global_id( 0 );

	point  p = dPobj[gid];
	vector v = dVel[gid];

	point  pp = p + v*DT + G*(point)( .5*DT*DT );
	vector vp = v + G*DT;
	pp.w = 1.;
	vp.w = 0.;


	// Check not inside sphere 1 (Outer Sphere)
	if( !IsInsideSphere( pp, Sphere1 ) )
	{
		vp = BounceSphere( p, v, Sphere1 );
		pp = p + vp*DT + G*(point)( .5*DT*DT );

		dCobj[gid].r = 0.89;
		dCobj[gid].g = 0.47;
		dCobj[gid].b = 0.20;
		dCobj[gid].a = 1.0;		// Color orange on bounce

	}

	
	// Sphere if inside sphere 2	(Inner Sphere)
	if( IsInsideSphere( pp, Sphere2 ) )
	{
		vp = BounceSphere( p, v, Sphere2 );
		pp = p + vp*DT + G*(point)( .5*DT*DT );	
		

		dCobj[gid].r = 0.90;
		dCobj[gid].g = 0.91;
		dCobj[gid].b = 0.98;
		dCobj[gid].a = 0.1;		// Color silver on bounce

	}

	dPobj[gid] = pp;
	dVel[gid]  = vp;
}




