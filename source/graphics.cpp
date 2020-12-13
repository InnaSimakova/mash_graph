#include "stdafx.h"    

using namespace std;

double wd = 800, ht = 800;
bool shadow_picture_flag = false, real_image_flag = false;
int effects[] = {1, 1, 1, 1, 1};

GLfloat ambient[] = {0.7f, 0.7f, 0.7f}; 
GLfloat diffuse[] = {1.0f, 1.0f, 1.0f}; 
GLfloat position[] = {0.0f, 4.0f, 0.0f, 1.0f};
GLfloat spot[] = {0.0f, 1.0f, 0.0f};
 
GLuint texture[2], nurbs_list;   
GLUquadricObj *q = gluNewQuadric();

struct camera
{
	GLfloat x, y, z;
	GLfloat zoom, rot;
	GLfloat dr;
	GLfloat dy;
	GLfloat dz;

	camera()
    { 
		rot = -45.0;
		dr = 0.5f;
		dy = 0.1f;
		dz = 0.2f;
		zoom = -3.0f;
        x = cos(rot/180.0 * M_PI) * zoom;
        y = 1.0;
        z = sin(rot/180.0 * M_PI) * zoom;
    }
};

struct camera cam;

struct part_type{

	GLfloat x, y, z;
	GLfloat dx, dy, dz;
	GLfloat speed, size;

	bool operator<(const part_type& e)const;
};

bool part_type::operator<(const part_type& e)const{

	return 
		fabs(x - cam.x) + fabs(y - cam.y) + fabs(z - cam.z) > 
		fabs(e.x - cam.x) + fabs(e.y - cam.y) + fabs(e.z - cam.z);

}

struct poly_type
{
	Vector3D p1, p2, p3, p4;

	poly_type(Vector3D s1, Vector3D s2, Vector3D s3, Vector3D s4){ p1 = s1; p2 = s2; p3 = s3; p4 = s4; }

	bool operator<(const poly_type& e)const;
};

bool poly_type::operator<(const poly_type& e)const{

	return 
		fabs((p1.x + p2.x + p3.x + p4.x)/4.0 - cam.x) + fabs((p1.y + p2.y + p3.y + p4.y)/4.0 - cam.y) + fabs((p1.z + p2.z + p3.z + p4.z)/4.0 - cam.z) > 
		fabs((e.p1.x + e.p2.x + e.p3.x + e.p4.x)/4.0 - cam.x) + fabs((e.p1.y + e.p2.y + e.p3.y + e.p4.y)/4.0 - cam.y) + fabs((e.p1.z + e.p2.z + e.p3.z + e.p4.z)/4.0 - cam.z);

}

vector<part_type> parts; // билборды
vector<poly_type> parts_obj; // полупрозрачный объект

void draw_base()
{
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture[0]);

	GLfloat mat_diffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat mat_specular[] = {0.8f, 0.8f, 0.8f, 1.0f};
    GLfloat mat_shininess[] = {0.0f};

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	glPushMatrix();
	glBegin(GL_QUADS);

	glNormal3f(0.0, 1.0, 0.0);
	glTexCoord2f(0.0f, 5.0f);  
	glVertex3f(-2.0, 0.0, 2.0);
	glTexCoord2f(0.0f, 0.0f);  
	glVertex3f(-2.0, 0.0,-2.0);
	glTexCoord2f(5.0f, 0.0f); 
	glVertex3f( 2.0, 0.0,-2.0);
	glTexCoord2f(5.0f, 5.0f);  
	glVertex3f( 2.0, 0.0, 2.0);
	
	glEnd();                     
	glPopMatrix();

	glDisable(GL_TEXTURE_2D);
}

void draw_mirror()
{
	GLUquadric *qobj = gluNewQuadric();
	glPushMatrix();

	glRotatef(90.0, 0.0, 1.0, 0.0);
	glTranslatef(0.0, 1.0, 0.0);
	gluDisk(qobj, 0.0, 1.0, 20, 20);

	glPopMatrix();
}

void sub(void (*draw)())
{
	glClear ( GL_STENCIL_BUFFER_BIT );
	
	glEnable      ( GL_CULL_FACE );
	glEnable      ( GL_DEPTH_TEST );
	glColorMask   ( GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE );
	glDepthMask   ( GL_FALSE );
	glCullFace    ( GL_BACK );
	glStencilOp   ( GL_KEEP, GL_KEEP, GL_INCR );
	glStencilFunc ( GL_ALWAYS, 0, 0 );
	glEnable      ( GL_STENCIL_TEST );

	draw ();

	glCullFace    ( GL_FRONT );
	glDepthMask   ( GL_TRUE  );
	glStencilOp   ( GL_KEEP, GL_KEEP, GL_KEEP );
	glStencilFunc ( GL_NOTEQUAL, 0, 1 );
	glDepthFunc   ( GL_GREATER );
	
	draw ();
	
	glDisable     ( GL_STENCIL_TEST );
	glDepthFunc   ( GL_LESS );
	glCullFace    ( GL_BACK );
}

void clip(void (*draw)())
{
	glCullFace    ( GL_FRONT );
	
	draw ();
	
	glCullFace    ( GL_BACK );
}

void draw_1()
{
	glPushMatrix();

	glutSolidSphere(0.5, 20, 20);
	glPopMatrix();
}

void draw_2()
{
	glPushMatrix();

	glTranslatef(0.0f, 0.5 ,0.0f);
	glutSolidSphere(0.4, 20, 20);
		
	glPopMatrix();
}

void draw_3()
{
	glPushMatrix();

	glTranslatef(0.0f, 0.2f ,0.0f);
	glutSolidCube(0.3);

	glPopMatrix();
}

void draw_csg()
{
	glPushMatrix();
	glTranslatef(1.0, 1.0, 0.0);
	glDisable(GL_TEXTURE_2D);

	GLfloat mat_diffuse[] = {0.6f, 0.6f, 0.6f, 1.0f};
	GLfloat mat_specular[] = {0.4f, 0.4f, 0.4f, 1.0f};
	GLfloat mat_shininess[] = {50.0f};

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glColor4f(0.2, 0.6, 0.8, 1.0);

	glPushMatrix();

	glEnable     ( GL_CULL_FACE );
	glCullFace   ( GL_BACK );
	glDepthFunc  ( GL_LESS );
	
	glColorMask   ( 0,0,0,0 );		
	
	draw_1 ();

	sub   ( draw_2 );
	sub   ( draw_3 );		
	clip  ( draw_1 );
	
	glColorMask ( 1,1,1,1 );
	glDepthFunc ( GL_EQUAL );

	glCullFace ( GL_BACK );

	draw_1 ();
	
	glCullFace ( GL_FRONT );

	glPushMatrix();
	glScalef(-1, -1, -1);
	glLightfv(GL_LIGHT0, GL_POSITION, position);      
	glPopMatrix();
	
	draw_2 ();
	draw_3 ();

	glDepthFunc ( GL_LEQUAL );

	glDisable(GL_STENCIL_TEST);
	glDisable(GL_CULL_FACE);

	glPopMatrix();

	glEnable(GL_TEXTURE_2D);
	glPopMatrix();
	glLightfv(GL_LIGHT0, GL_POSITION, position);      

}

void init_parts ()
{
	srand(time(NULL));

	for(int i = 0; i < 100; i++) 
	{
		struct part_type part;
		part.size = rand()/(float)RAND_MAX * 0.05;
		part.speed = 0.005 + rand()/(float)RAND_MAX*0.005;
		part.x = -2 + rand()/(float)RAND_MAX * 4.0;
		part.y = rand()/(float)RAND_MAX * 4.0;
		part.z = -2 + rand()/(float)RAND_MAX * 4.0;
		part.dx = rand() % 10;
		part.dy = rand() % 10;
		part.dz = rand() % 10;

		parts.push_back(part);
	}

}

void billboardSphericalBegin( float camX, float camY, float camZ,
                              float objPosX, float objPosY, float objPosZ)
{
  float lookAt[3],objToCamProj[3],objToCam[3],upAux[3];
  float angleCosine;

  glPushMatrix();

  objToCamProj[0] = camX - objPosX ;
  objToCamProj[1] = 0.0;
  objToCamProj[2] = camZ - objPosZ ;

  lookAt[0] = 0.0;
  lookAt[1] = 0.0;
  lookAt[2] = 1.0;

  mathsNormalize(objToCamProj);
  mathsCrossProduct(upAux,lookAt,objToCamProj);
  angleCosine = mathsInnerProduct(lookAt,objToCamProj);

  if ((angleCosine <= 1.0) && (angleCosine >= -1.0))
    glRotatef(acos(angleCosine)*180.0/M_PI,upAux[0], upAux[1], upAux[2]);

  objToCam[0] = camX - objPosX;
  objToCam[1] = camY - objPosY;
  objToCam[2] = camZ - objPosZ;

  mathsNormalize(objToCam);

  angleCosine = mathsInnerProduct(objToCamProj,objToCam);

  if ((angleCosine <= 1.0) && (angleCosine >= -1.0))
    if (objToCam[1] < 0.0)
      glRotatef(acos(angleCosine)*180.0/M_PI, 1.0, 0.0, 0.0);
    else
      glRotatef(acos(angleCosine)*180.0/M_PI, -1.0, 0.0, 0.0);
}

void billboardEnd()
{
  glPopMatrix();
}

void draw_parts()
{
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture[0]);
	if(!shadow_picture_flag)
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_BLEND);

		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		GLfloat mat_diffuse[] =	{0.9f, 0.9f, 0.9f, 0.5f};
		GLfloat mat_specular[] = {1.0f, 1.0f, 1.0f, 0.5f};
		GLfloat mat_shininess[] = {150.0f};

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
		glColor4f(1,1,1,0.5);
	}

	sort(parts.begin(), parts.end());

	for (unsigned int i = 0; i < parts.size(); i++)
	{
		glPushMatrix();
		glTranslatef(parts[i].x, parts[i].y, parts[i].z);

		
		if(!shadow_picture_flag) billboardSphericalBegin(cam.x, cam.y, cam.z,
                              parts[i].x, parts[i].y, parts[i].z);

		else
		{
			glPushMatrix();
			glRotatef(90.0, 1.0, 0.0, 0.0);
		}

		glPushMatrix();
		gluDisk(q, 0.0, parts[i].size, 6, 6);
		glPopMatrix();

		if(!shadow_picture_flag) billboardEnd();
		else glPopMatrix();

		glPopMatrix();

	}

	if(!shadow_picture_flag)
	{
		glColor4f(1,1,1,1);
		glDisable(GL_BLEND);
		glEnable(GL_LIGHTING);
	}
}

void parts_next_frame()
{
	for (unsigned int i = 0; i < parts.size(); i++)
	{
		parts[i].y -= parts[i].speed;

		if (parts[i].y < 0)
		{
				parts[i].y = 4;
				parts[i].x = -2 + rand()/(float)RAND_MAX * 4.0;
				parts[i].z =-2 + rand()/(float)RAND_MAX * 4.0;
				parts[i].size = rand()/(float)RAND_MAX * 0.05;
				parts[i].speed = 0.005 + rand()/(float)RAND_MAX*0.005;
		}
	}
}


Vector3D figure_parametrization(float u, float v)
{
   Vector3D p;

   p.x = (2 + u * cos(v)) * sin(2 * M_PI * u);
   p.y = (2 + u * cos(v)) * cos(2 * M_PI * u) + 2 * u;
   p.z = u * sin(v);

   p /= 3.0;
   p += Vector3D(-1.0, 0.5, -1.5);

   return p;
}

void init_parts_obj()
{
	parts_obj.clear();

	for(int i = 0; i < 50; i++)
	{
		for(int j = 0; j < 50; j++)
		{
			float u  = i / 50.0;
			float u_next = (i + 1) / 50.0;
			float v  = j * M_PI * 2.0 / 50.0;
			float v_next = (j + 1) * M_PI * 2.0 / 50.0;

			Vector3D p1 = figure_parametrization(u, v);
			Vector3D p2 = figure_parametrization(u_next, v);
			Vector3D p3 = figure_parametrization(u_next, v_next);
			Vector3D p4 = figure_parametrization(u, v_next);

			struct poly_type q(p1, p2, p3, p4);

			parts_obj.push_back(q);
      }
    }
}

void draw_parts_obj()
{
	if(!shadow_picture_flag)
	{
		GLfloat mat_diffuse[] = {0.6f, 0.6f, 0.6f, 1.0f};
		GLfloat mat_specular[] = {0.4f, 0.4f, 0.4f, 1.0f};
		GLfloat mat_shininess[] = {50.0f};

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

		glDisable(GL_TEXTURE_2D);
		glEnable(GL_BLEND);

		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	}

	sort(parts_obj.begin(), parts_obj.end());
	glEnable(GL_NORMALIZE);
	glPushMatrix();

	for(unsigned int i = 0; i < parts_obj.size(); i++)
	{
		Vector3D n = 
			(Vector3D(parts_obj[i].p2.x, parts_obj[i].p2.y, parts_obj[i].p2.z) - 
			Vector3D(parts_obj[i].p3.x, parts_obj[i].p3.y, parts_obj[i].p3.z)) ^ 
			(Vector3D(parts_obj[i].p1.x, parts_obj[i].p1.y, parts_obj[i].p1.z) - 
			Vector3D(parts_obj[i].p2.x, parts_obj[i].p2.y, parts_obj[i].p2.z)); 
			
		n.normalize();
		glNormal3f(n.x, n.y, n.z);

		glBegin(GL_QUADS);
	 
		if(!shadow_picture_flag) glColor4f(parts_obj[i].p1.x + 1.0, parts_obj[i].p1.y - 0.5, parts_obj[i].p1.z + 1.5, 0.6f);
		glVertex3f(parts_obj[i].p1.x, parts_obj[i].p1.y, parts_obj[i].p1.z);
		if(!shadow_picture_flag) glColor4f(parts_obj[i].p2.x + 1.0, parts_obj[i].p2.y - 0.5, parts_obj[i].p2.z + 1.5, 0.6f);
		glVertex3f(parts_obj[i].p2.x, parts_obj[i].p2.y, parts_obj[i].p2.z);
		if(!shadow_picture_flag) glColor4f(parts_obj[i].p3.x + 1.0, parts_obj[i].p3.y - 0.5, parts_obj[i].p3.z + 1.5, 0.6f);
		glVertex3f(parts_obj[i].p3.x, parts_obj[i].p3.y, parts_obj[i].p3.z);
		if(!shadow_picture_flag) glColor4f(parts_obj[i].p4.x + 1.0, parts_obj[i].p4.y - 0.5, parts_obj[i].p4.z + 1.5, 0.6f);
		glVertex3f(parts_obj[i].p4.x, parts_obj[i].p4.y, parts_obj[i].p4.z);

		glEnd();
	}

	glPopMatrix();
	glDisable(GL_NORMALIZE);

	if(!shadow_picture_flag)
	{
		glColor4f(1.0, 1.0, 1.0, 1.0);
		glEnable(GL_TEXTURE_2D);
		glDisable(GL_BLEND);
	}
}


float spline_param = 0.0;

void init_spline_obj()
{
	GLUnurbsObj *nurb = gluNewNurbsRenderer(); 
	GLfloat control_points[4][4][3] =
	{
	{
		{0.0f, 0.0f, 0.0f},
		{0.0f, 0.0f, 0.2f},
		{0.0f, 0.0f, 0.4f},
		{0.0f, 0.0f, 0.6f},
	},
	{
		{- cos(spline_param) * 1.5, 1.5f, 0.0f},
		{- cos(spline_param) * 1.5, 1.5f, 0.2f},
		{- cos(spline_param) * 1.5, 1.5f, 0.4f},
		{- cos(spline_param) * 1.5, 1.5f, 0.6f},
	},
	{
		{cos(spline_param) * 1.5, 2.5f, 0.0f},
		{cos(spline_param) * 1.5, 2.5f, 0.2f},
		{cos(spline_param) * 1.5, 2.5f, 0.4f},
		{cos(spline_param) * 1.5, 2.5f, 0.6f},
	},
	{
		{0.0f, 4.0f, 0.0f},
		{0.0f, 4.0f, 0.2f},
		{0.0f, 4.0f, 0.4f},
		{0.0f, 4.0f, 0.6f},
	}};

    gluNurbsProperty(nurb, GLU_SAMPLING_TOLERANCE, 25.0);
    gluNurbsProperty(nurb, GLU_DISPLAY_MODE, GLU_FILL);

	nurbs_list = glGenLists(1);
	glNewList(nurbs_list, GL_COMPILE);

	GLfloat knots[8] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
	GLfloat texture_points[4][4][3] =
	{
	{
		{0.0f, 1.0f, 0.0f},
		{0.0f, 0.66f, 0.0f},
		{0.0f, 0.33f, 0.0f},
		{0.0f, 0.0f, 0.0f},
	},
	{
		{0.2f, 1.0f, 0.0f},
		{0.2f, 0.66f, 0.0f},
		{0.2f, 0.33f, 0.0f},
		{0.2f, 0.0f, 0.0f},
	},
	{
		{0.8f, 1.0f, 0.0f},
		{0.8f, 0.66f, 0.0f},
		{0.8f, 0.33f, 0.0f},
		{0.8f, 0.0f, 0.0f},
	},
	{
		{1.0f, 1.0f, 0.0f},
		{1.0f, 0.66f, 0.0f},
		{1.0f, 0.33f, 0.0f},
		{1.0f, 0.0f, 0.0f},
	}};

	glPushMatrix();
	glEnable(GL_AUTO_NORMAL);

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture[1]);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	
	glPushMatrix();

	glTranslatef(0.0, 3.0, -2.0);
	glRotatef(90.0, 1.0, 0.0, 0.0);
	glNormal3f(1.0, 0.0, 0.0);

    gluBeginSurface(nurb);
    gluNurbsSurface(nurb, 8, knots, 8, knots, 4 * 3, 3, &texture_points[0][0][0], 4, 4, GL_MAP2_TEXTURE_COORD_3);
    gluNurbsSurface(nurb, 8, knots, 8, knots, 4 * 3, 3, &control_points[0][0][0], 4, 4, GL_MAP2_VERTEX_3);

    gluEndSurface(nurb);
	glDisable(GL_TEXTURE_2D);

    glPopMatrix();
	glPopMatrix();

	glEndList();
}

void draw_spline()
{
	init_spline_obj();
	glPushMatrix();
	glCallList(nurbs_list);
	glPopMatrix();
}

void draw_teapot()
{
	if(!shadow_picture_flag) glColor4f(1.0, 1.0, 1.0, 1.0);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture[0]);
	glPushMatrix();
	glTranslatef(-1.0, 1.0, 0.0);
	glutSolidTeapot(0.3);
	glPopMatrix();
	glDisable(GL_TEXTURE_2D);
}

int load_textures()
{
	AUX_RGBImageRec *texture_image[2];        

	texture_image[0] = auxDIBImageLoad(_T("../data/base.bmp"));
	texture_image[1] = auxDIBImageLoad(_T("../data/flag.bmp"));

	if(texture_image[0] == NULL) return false;
	if(texture_image[1] == NULL) return false;
	glGenTextures(2, &texture[0]);

	glBindTexture(GL_TEXTURE_2D, texture[0]);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_NEAREST);
	gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texture_image[0]->sizeX, texture_image[0]->sizeY, GL_RGB, GL_UNSIGNED_BYTE, texture_image[0]->data); 

	glBindTexture(GL_TEXTURE_2D, texture[1]);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, texture_image[1]->sizeX, texture_image[1]->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, texture_image[1]->data);

	return true;
}

int init_gl()       
{
	init_parts();
	init_spline_obj();
	init_parts_obj();
	load_textures(); 

	GLfloat mat_diffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat mat_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat mat_shininess[] = {50.0f};

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	glShadeModel(GL_SMOOTH);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glColor4f(1.0, 1.0, 1.0, 1.0);
	glClearDepth(1.0f);   
	glClearStencil(0);  
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL); 

	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
 	glLightfv(GL_LIGHT0, GL_POSITION, position);      
		
	glEnable(GL_LIGHT0);            
	glEnable(GL_LIGHTING); 

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    gluQuadricNormals(q, GL_SMOOTH);
    gluQuadricTexture(q, GL_TRUE);
    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
    glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);

	return true; 
}

void shadowmatrix(GLfloat matrix[4][4], GLfloat plane[4], GLfloat position[4])
{
	GLfloat dot;

	dot = plane[0] * position[0] +
		plane[1] * position[1] +
		plane[2] * position[2] +
		plane[3];

	matrix[0][0] = dot - position[0] * plane[0];
	matrix[1][0] = 0.f - position[0] * plane[1];
	matrix[2][0] = 0.f - position[0] * plane[2];
	matrix[3][0] = 0.f - position[0] * plane[3];

	matrix[0][1] = 0.f - position[1] * plane[0];
	matrix[1][1] = dot - position[1] * plane[1];
	matrix[2][1] = 0.f - position[1] * plane[2];
	matrix[3][1] = 0.f - position[1] * plane[3];

	matrix[0][2] = 0.f - position[2] * plane[0];
	matrix[1][2] = 0.f - position[2] * plane[1];
	matrix[2][2] = dot - position[2] * plane[2];
	matrix[3][2] = 0.f - position[2] * plane[3];

	matrix[0][3] = 0.f - position[3] * plane[0];
	matrix[1][3] = 0.f - position[3] * plane[1];
	matrix[2][3] = 0.f - position[3] * plane[2];
	matrix[3][3] = dot - position[3] * plane[3];
}

void draw_shadow()
{
    GLfloat matrix[4][4], plane[] = {0.0f, 1.0f, 0.0f, 0.0f};

	shadow_picture_flag = true;

    glEnable(GL_STENCIL_TEST);

	shadowmatrix(matrix, plane, position);

	glClear(GL_STENCIL_BUFFER_BIT);
	glClearStencil(0);
	glStencilFunc(GL_ALWAYS, 0x1, 0xffffffff);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glColorMask(0, 0, 0, 0);

	draw_base();
	
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);

	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glColorMask(1, 1, 1, 1);
	glColor4f(0.0f, 0.0f, 0.0f, 0.2f);

	glStencilOp(GL_KEEP,GL_DECR,GL_DECR);
	glStencilFunc(GL_EQUAL,1,1);

	glPushMatrix();

	glMultMatrixf((GLfloat *) matrix);
	if(effects[1]) draw_spline();
	if(effects[2]) draw_parts_obj();
	if(effects[3]) draw_parts();
	draw_teapot();

	if(effects[4])
	{
		glPushMatrix();
		glTranslatef(1.0, 1.0, 0.0);
		glutSolidSphere(0.5, 20, 20);
		glPopMatrix();
	}
	glPopMatrix();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
   
	glDisable(GL_BLEND);
	glDisable(GL_STENCIL_TEST);

   shadow_picture_flag = false;
}


int draw_scene()
{
	double eqr[] = {1.0f, 0.0f, 0.0f, -0.001f}; 

	if(effects[0] && cam.x <= 0)
	{

		glClear(GL_STENCIL_BUFFER_BIT);
		glEnable(GL_STENCIL_TEST);
		glClearStencil(0);
		glStencilFunc(GL_ALWAYS, 1, 0xffffffff);
		glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
		glColorMask(0,0,0,0);
		glDisable(GL_DEPTH_TEST);
	
		draw_mirror();

		glEnable(GL_DEPTH_TEST);
		glColorMask(1, 1, 1, 1);
		glStencilFunc(GL_NOTEQUAL, 0, 0xffffffff);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

		glPushMatrix();

		glEnable(GL_CLIP_PLANE0);
		glClipPlane(GL_CLIP_PLANE0, eqr);
		glScalef(-1.0f, 1.0f, 1.0f);
		glLightfv(GL_LIGHT0, GL_POSITION, position); 

		cam.x = -cam.x;

		draw_base();
		draw_teapot();
	
		if(effects[1]) draw_spline();
		if(effects[2]) draw_parts_obj();
		if(effects[3]) draw_parts ();

		glPopMatrix();   
		cam.x = -cam.x;
		glDisable(GL_CLIP_PLANE0);

		glLightfv(GL_LIGHT0, GL_POSITION, position); 

		glClear(GL_STENCIL_BUFFER_BIT);
		glEnable(GL_STENCIL_TEST);
		glClearStencil(0);
		glStencilFunc(GL_ALWAYS, 1, 0xffffffff);
		glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
		glColorMask(0,0,0,0);

		draw_mirror();

		glColorMask(1,1,1,1);
		glDisable(GL_STENCIL_TEST);

	}
	real_image_flag = true;

	draw_base();
	draw_shadow();
	draw_teapot();
	if(effects[1]) draw_spline();
	if(effects[2]) draw_parts_obj();
	if(effects[3]) draw_parts ();
	if(effects[0] && cam.x > 0)
	{
		glPushMatrix();
		glDisable(GL_LIGHTING);
		glColor3f(0.0, 0.0, 0.0);
		draw_mirror();
		glEnable(GL_LIGHTING);
		glColor3f(1.0, 1.0, 1.0);
		glPopMatrix();
	}
	if(effects[4]) draw_csg();

	real_image_flag = false;

	return true; 
}

int FPS = 0;

void calculate_frame_rate()
{
    static float framesPerSecond = 0.0f; 
    static float lastTime = 0.0f;         

	float currentTime = glutGet(GLUT_ELAPSED_TIME) * 0.001f;

    ++framesPerSecond;

    if(currentTime - lastTime > 1.0f)
    {
        lastTime = currentTime;
        FPS = framesPerSecond;
        framesPerSecond = 0;
    }
}


void display()
{
    glMatrixMode( GL_MODELVIEW ); 
	glLoadIdentity();

	gluLookAt(
        cam.x, cam.y, cam.z,
        0.0, 1.0, 0.0,
        0.0f, 1.0f, 0.0f); 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	glLightfv(GL_LIGHT0, GL_POSITION, position);

	draw_scene();
    calculate_frame_rate();
	char* buff = (char *)malloc(10 * sizeof(char));
	sprintf(buff, "FPS: %d", FPS);
	glutSetWindowTitle(buff);

	glutSwapBuffers();   
}

void reshape(int w, int h)
{
	if(h == 0) h = 1;

	glViewport(0, 0, w, h); 
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(70.0f, (GLfloat)w/(GLfloat)h, 0.1f, 100.0f );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();              
}

void process_keyboard(unsigned char key, int x, int y)  
{
	if(key == '\033')
      std::exit(0);

	if(key == '1')
		effects[0] = 1 - effects[0];
	if(key == '2')
		effects[1] = 1 - effects[1];
	if(key == '3')
		effects[2] = 1 - effects[2];
	if(key == '4')
		effects[3] = 1 - effects[3];
	if(key == '5')
		effects[4] = 1 - effects[4];
	if(key == '6')
		effects[5] = 1 - effects[5];
	if(key == '7')
		effects[6] = 1 - effects[6];
}

void process_spec(int key, int x, int y)
{
	if(key == GLUT_KEY_DOWN)
	{
			cam.y -= cam.dy; 
    }    
	if(key == GLUT_KEY_UP)
	{
			cam.y += cam.dy;
    }
	if(key == GLUT_KEY_RIGHT)
	{
		cam.rot -= cam.dr;
		cam.x = cos(cam.rot/180.0 * M_PI) * cam.zoom;
        cam.z = sin(cam.rot/180.0 * M_PI) * cam.zoom;
	}
	if(key == GLUT_KEY_LEFT)
	{
		cam.rot += cam.dr;
		cam.x = cos(cam.rot/180.0 * M_PI) * cam.zoom;
        cam.z = sin(cam.rot/180.0 * M_PI) * cam.zoom;
    }
	if(key == GLUT_KEY_PAGE_UP)
	{
		cam.zoom += cam.dz;
		cam.x = cos(cam.rot/180.0 * M_PI) * cam.zoom;
        cam.z = sin(cam.rot/180.0 * M_PI) * cam.zoom;
    }
	if(key == GLUT_KEY_PAGE_DOWN)
	{
		cam.zoom -= cam.dz;
		cam.x = cos(cam.rot/180.0 * M_PI) * cam.zoom;
        cam.z = sin(cam.rot/180.0 * M_PI) * cam.zoom;
    }

	glutPostRedisplay();
}

void Timer(int value)
{
	parts_next_frame();
	spline_param += 0.01;
	glutTimerFunc(20, Timer, 0);
}

void idle()
{
	glutPostRedisplay();

}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
    glutInitWindowSize(wd, ht);
	glutCreateWindow(argv[0]);
    glutReshapeFunc(reshape);
	glutKeyboardFunc(process_keyboard);
	glutSpecialFunc(process_spec);
    glutDisplayFunc(display);
	glutTimerFunc(20, Timer, 0);
	glutIdleFunc(idle);
    init_gl();
    glutMainLoop();

    return 0;             
}