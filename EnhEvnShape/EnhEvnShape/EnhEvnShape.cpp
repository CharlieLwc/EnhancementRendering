// EnhEvnShape.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

/*
Szymon Rusinkiewicz
Princeton University

mesh_view.cc
Simple viewer
*/
#define GLUT_DISABLE_ATEXIT_HACK
#include "../TriMesh2/TriMesh.h"
#include "../TriMesh2/XForm.h"
#include "../TriMesh2/GLCamera.h"
#include "../TriMesh2/ICP.h"
#include "../TriMesh2/strutil.h"


#include "barItems.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <string.h>
#include"../GL/glut.h"


TwBar *bar;
vector<TriMesh *> meshes;




vector<char*> picName;



SHLighting shLight(meshes, "grace_cross", 2);
Transport transport(meshes, "building_cross", 8, 0.0075f);
RadianceScaling radianceScaling;


// Globals


vector<xform> xforms;
vector<bool> visible;
vector<string> filenames;

TriMesh::BSphere global_bsph;
xform global_xf;
GLCamera camera;





int mouseXInc;
int mouseYInc;

bool grab_only = false;


int current_mesh = -1;


// Signal a redraw
void need_redraw()
{ 
	glutPostRedisplay();
}


// Clear the screen
void cls()
{
	glDisable(GL_DITHER);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_NORMALIZE);
	glDisable(GL_LIGHTING);
	glDisable(GL_NORMALIZE);
	glDisable(GL_COLOR_MATERIAL);
	if (ISdraw_white_bg())
		glClearColor(1, 1, 1, 0);
	else
		glClearColor(0.08f, 0.08f, 0.08f, 0);
	glClearDepth(1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}


// Set up lights and materials
void setup_lighting(int id)
{
	Color c(1.0f);
	if (ISdraw_falsecolor())
		c = Color::hsv(-3.88f * id, 0.6f + 0.2f * sin(0.42f * id), 1);
	

	glColor3fv(c);

	if (ISdisableLighting()) {
		glDisable(GL_LIGHTING);
		return;
	}


	GLfloat mat_specular[4] = { 0.18f, 0.18f, 0.18f, 0.18f };
	if (!ISdraw_shiny()) {
		mat_specular[0] = mat_specular[1] =
			mat_specular[2] = mat_specular[3] = 0.0f;
	}
	GLfloat mat_shininess[] = { 64 };
	GLfloat global_ambient[] = { 0.02f, 0.02f, 0.05f, 0.05f };
	GLfloat light0_ambient[] = { 0, 0, 0, 0 };
	GLfloat light0_diffuse[] = { 0.85f, 0.85f, 0.8f, 0.85f };
	if (current_mesh >= 0 && id != current_mesh) {
		light0_diffuse[0] *= 0.5f;
		light0_diffuse[1] *= 0.5f;
		light0_diffuse[2] *= 0.5f;
	}
	GLfloat light1_diffuse[] = { -0.01f, -0.01f, -0.03f, -0.03f };
	GLfloat light0_specular[] = { 0.85f, 0.85f, 0.85f, 0.85f };
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, ISdraw_2side());
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);


}

// Draw triangle strips.  They are stored as length followed by values.
void draw_tstrips(const TriMesh *themesh, bool useFaceNormal)
{

	drawEyedirOnScreen(global_bsph, global_xf);

	static bool use_glArrayElement = false;
	static bool tested_renderer = false;
	if (!tested_renderer) {
		use_glArrayElement = !!strstr(
			(const char *) glGetString(GL_RENDERER), "Intel");
		tested_renderer = true;
	}

	const int *t = &themesh->tstrips[0];
	const int *end = t + themesh->tstrips.size();


	if (useFaceNormal)
	{

		int faceIndex = 0;
		vec3 temp(1.0, 0.0, 0.0);
		while (likely(t < end)) {


			glBegin(GL_TRIANGLES);

			int striplen = *t++;


			int i = 3;
			for (i = 3; i < striplen; i+=2)
			{
			
				glNormal3fv(themesh->faceNormals[faceIndex]);
				glVertex3fv(themesh->vertices[*t++]);
				glNormal3fv(themesh->faceNormals[faceIndex]);
				glVertex3fv(themesh->vertices[*t++]);
				glNormal3fv(themesh->faceNormals[faceIndex]);
				glVertex3fv(themesh->vertices[*t]);

				faceIndex++;
				glNormal3fv(themesh->faceNormals[faceIndex]);
				glVertex3fv(themesh->vertices[*t--]);
				glNormal3fv(themesh->faceNormals[faceIndex]);
				glVertex3fv(themesh->vertices[*t++]);
				t++;
				glNormal3fv(themesh->faceNormals[faceIndex]);
				glVertex3fv(themesh->vertices[*t--]);
				faceIndex++;

				/*
				glArrayElement(*t++);
				glArrayElement(*t++);
				glArrayElement(*t);

				glArrayElement(*t--);
				glArrayElement(*t++);
				t++;
				glArrayElement(*t--);
				*/
			}

			if (i == striplen)
			{
				glNormal3fv(themesh->faceNormals[faceIndex]);
				glVertex3fv(themesh->vertices[*t++]);
				glNormal3fv(themesh->faceNormals[faceIndex]);
				glVertex3fv(themesh->vertices[*t++]);
				glNormal3fv(themesh->faceNormals[faceIndex]);
				glVertex3fv(themesh->vertices[*t++]);
				faceIndex++;
			}
			else
			{
				t++;
				t++;
			}

			glEnd();
			

		}
	}
	else
	{
		if (use_glArrayElement) {
			while (likely(t < end)) {
				glBegin(GL_TRIANGLE_STRIP);
				int striplen = *t++;
				for (int i = 0; i < striplen; i++)
					glArrayElement(*t++);
				glEnd();
			}
		}
		else {
			while (likely(t < end)) {
				int striplen = *t++;
				glDrawElements(GL_TRIANGLE_STRIP, striplen, GL_UNSIGNED_INT, t);
				t += striplen;
			}
		}
	}



}


// Draw the mesh
void draw_mesh(int meshIndex)
{
	TriMesh *themesh = meshes[meshIndex];



	glPushMatrix();

	glMultMatrixd(xforms[meshIndex]);

	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);

	if (ISdraw_2side()) {
		glDisable(GL_CULL_FACE);
	} else {
		glCullFace(GL_BACK);
		glEnable(GL_CULL_FACE);
	}

	// Vertices
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT,
		sizeof(themesh->vertices[0]),
		&themesh->vertices[0][0]);







	// Colors
	if (!themesh->colors.empty() && useColor()) {
		glEnableClientState(GL_COLOR_ARRAY);
		glColorPointer(3, GL_FLOAT,
			sizeof(themesh->colors[0]),
			&themesh->colors[0][0]);
	}
	else {
		glDisableClientState(GL_COLOR_ARRAY);
	}

	// Main drawing pass
	if (ISdraw_points() || themesh->tstrips.empty()) {
		// No triangles - draw as points
		glPointSize(float(getdraw_point_size()));
		glDrawArrays(GL_POINTS, 0, themesh->vertices.size());
		glPopMatrix();
		return;
	}

	if (ISdraw_edges() || ISdraw_ridge()) {
		glPolygonOffset(10.0f, 10.0f);
		glEnable(GL_POLYGON_OFFSET_FILL);
	}


	meshSegPreparetion();
	normalPreparetion();
	curvPreparetion();
	surfaceReliefPreparetion(global_bsph, global_xf);
	curveAndPlanePreparetion(global_bsph, global_xf);
	PRTPreparetion(global_bsph);

	void(*func)(const TriMesh *, const bool);
	func = draw_tstrips;
	if (!GLSLPreparetion(global_bsph, camera, func, global_xf))
	{
		if (!DrawLabSpace(global_bsph))
			draw_tstrips(themesh, IsFaceNormal());




	}
	
	


	glDisable(GL_POLYGON_OFFSET_FILL);

	// Edge drawing pass
	if (ISdraw_edges()) {
		glPolygonMode(GL_FRONT, GL_LINE);
		glLineWidth(float(getdraw_line_width()));
		glDisableClientState(GL_COLOR_ARRAY);
		glDisable(GL_COLOR_MATERIAL);
		GLfloat global_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
		GLfloat light0_diffuse[] = { 0.8f, 0.8f, 0.8f, 0.0f };
		GLfloat light1_diffuse[] = { -0.2f, -0.2f, -0.2f, 0.0f };
		GLfloat light0_specular[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
		GLfloat mat_diffuse[4] = { 0.0f, 0.0f, 1.0f, 1.0f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
		glColor3f(0, 0, 1); // Used iff unlit
		draw_tstrips(themesh, IsFaceNormal());
		glPolygonMode(GL_FRONT, GL_FILL);
	}

	DrawRidgeLines();

	glPopMatrix();
}


// Draw the scene
void redraw()
{


	SLmatrix(global_xf);

	timestamp t = now();
	camera.setupGL(global_xf * global_bsph.center, global_bsph.r);

	glPushMatrix();
	glMultMatrixd(global_xf);
	cls();
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])
			continue;
		setup_lighting(i);
		draw_mesh(i);
	}

	glPopMatrix();

	TwDraw();
	glutSwapBuffers();



	if (grab_only) {
		void dump_image();
		dump_image();
		exit(0);
	}
	printf("\r                        \r%.1f msec.", 1000.0f * (now() - t));
	fflush(stdout);
}


// Update global bounding sphere.
void update_bsph()
{
	point boxmin(1e38f, 1e38f, 1e38f);
	point boxmax(-1e38f, -1e38f, -1e38f);
	bool some_vis = false;
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])	
			continue;
		some_vis = true;
		point c = xforms[i] * meshes[i]->bsphere.center;
		float r = meshes[i]->bsphere.r;
		for (int j = 0; j < 3; j++) {
			boxmin[j] = min(boxmin[j], c[j]-r);
			boxmax[j] = max(boxmax[j], c[j]+r);
		}
	}
	if (!some_vis)
		return;
	point &gc = global_bsph.center;
	float &gr = global_bsph.r;
	gc = 0.5f * (boxmin + boxmax);
	gr = 0.0f;
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])	
			continue;
		point c = xforms[i] * meshes[i]->bsphere.center;
		float r = meshes[i]->bsphere.r;
		gr = max(gr, dist(c, gc) + r);
	}
}


// Set the view...
void resetview()
{
	camera.stopspin();

	// Reload mesh xforms
	for (size_t i = 0; i < meshes.size(); i++)
		if (!xforms[i].read(xfname(filenames[i])))
			xforms[i] = xform();

	update_bsph();

	// Set camera to first ".camxf" if we have it...
	for (size_t i = 0; i < filenames.size(); i++) {
		if (global_xf.read(replace_ext(filenames[i], "camxf")))
			return;
	}

	// else default view
	global_xf = xform::trans(0, 0, -5.0f * global_bsph.r) *
		xform::trans(-global_bsph.center);
}


// Make some mesh current
void set_current(int i)
{
	camera.stopspin();
	if (false && i >= 0 && i < (int)meshes.size() && visible[i])
		current_mesh = i;
	else
		current_mesh = -1;
	need_redraw();
}


// Change visiblility of a mesh
void toggle_vis(int i)
{
	if (i >= 0 && i < (int)meshes.size())
		visible[i] = !visible[i];
	if (current_mesh == i && !visible[i])
		set_current(-1);
	update_bsph();
	need_redraw();
}


// Save the current image to a PPM file.
// Uses the next available filename matching filenamepattern
void dump_image()
{
	// Find first non-used filename
	const char filenamepattern[] = "img%d.ppm";
	int imgnum = 0;
	FILE *f;
	while (1) {
		char filename[1024];
		sprintf_s(filename, filenamepattern, imgnum++);
		fopen_s(&f,filename, "rb");
		if (!f) {
			fopen_s(&f,filename, "rb");
			printf("\n\nSaving image %s... ", filename);
			fflush(stdout);
			break;
		}
		fclose(f);
	}

	// Read pixels
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	GLint width = V[2], height = V[3];
	char *buf = new char[width*height*3];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(V[0], V[1], width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	// Flip top-to-bottom
	for (int i = 0; i < height/2; i++) {
		char *row1 = buf + 3 * width * i;
		char *row2 = buf + 3 * width * (height - 1 - i);
		for (int j = 0; j < 3 * width; j++)
			swap(row1[j], row2[j]);
	}

	// Write out file
	fprintf(f, "P6\n#\n%d %d\n255\n", width, height);
	fwrite(buf, width*height*3, 1, f);
	fclose(f);
	delete [] buf;

	printf("Done.\n\n");
}


// Save scan transforms
void save_xforms()
{
	for (size_t i = 0; i < xforms.size(); i++) {
		string xffile = xfname(filenames[i]);
		printf("Writing %s\n", xffile.c_str());
		xforms[i].write(xffile);
	}
}


// Save camera xform
void save_cam_xform()
{
	std::string camfile = replace_ext(filenames[0], "camxf");
	printf("Writing %s\n", camfile.c_str());
	global_xf.write(camfile);
}


// ICP
void do_icp(int n)
{
	camera.stopspin();

	if (current_mesh < 0 || current_mesh >= (int)meshes.size())
		return;
	if (n < 0 || n >= (int)meshes.size())
		return;
	if (!visible[n] || !visible[current_mesh] || n == current_mesh)
		return;
	ICP(meshes[n], meshes[current_mesh], xforms[n], xforms[current_mesh], 2);
	update_bsph();
	need_redraw();
}


// Handle mouse button and motion events
static unsigned buttonstate = 0;

void doubleclick(int button, int x, int y)
{
	// Render and read back ID reference image
	camera.setupGL(global_xf * global_bsph.center, global_bsph.r);
	glDisable(GL_BLEND);
	glDisable(GL_LIGHTING);
	glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glPushMatrix();
	glMultMatrixd(global_xf);
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])
			continue;
		glColor3ub((i >> 16) & 0xff,
			(i >> 8)  & 0xff,
			i        & 0xff);
		draw_mesh(i);
	}
	glPopMatrix();
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	y = int(V[1] + V[3]) - 1 - y;
	unsigned char pix[3];
	glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, pix);
	int n = (pix[0] << 16) + (pix[1] << 8) + pix[2];

	if (button == 0 || buttonstate == (1 << 0)) {
		// Double left click - select a mesh
		set_current(n);
	} else if (button == 2 || buttonstate == (1 << 2)) {
		// Double right click - ICP current to clicked-on
		do_icp(n);
	}
}

void mousehelperfunc(int x, int y)
{
	static const Mouse::button physical_to_logical_map[] = {
		Mouse::NONE, Mouse::ROTATE, Mouse::MOVEXY, Mouse::MOVEZ,
		Mouse::MOVEZ, Mouse::MOVEXY, Mouse::MOVEXY, Mouse::MOVEXY,
	};

	Mouse::button b = Mouse::NONE;
	if (buttonstate & (1 << 3))
		b = Mouse::WHEELUP;
	else if (buttonstate & (1 << 4))
		b = Mouse::WHEELDOWN;
	else if (buttonstate & (1 << 30))
		b = Mouse::LIGHT;
	else
		b = physical_to_logical_map[buttonstate & 7];

	if (current_mesh < 0) {
		camera.mouse(x, y, b,
			global_xf * global_bsph.center, global_bsph.r,
			global_xf, mouseXInc, mouseYInc);
	} else {
		xform tmp_xf = global_xf * xforms[current_mesh];
		camera.mouse(x, y, b,
			tmp_xf * meshes[current_mesh]->bsphere.center,
			meshes[current_mesh]->bsphere.r,
			tmp_xf, mouseXInc, mouseYInc);
		xforms[current_mesh] = inv(global_xf) * tmp_xf;
		update_bsph();
	}


	if (ISUsingSH() && (mouseXInc != 0 || mouseYInc != 0))
	{
		shLight.rote(mouseXInc, mouseYInc);
		shLight.setColor(meshes[0]->getColor());
	} 
}

bool mouseClickBar;

void mousemotionfunc(int x, int y)
{
	if(TwMouseMotion(x, y))
	{
		TwRefreshBar(bar);
		need_redraw();
		return;
	}

	mousehelperfunc(x,y);
	if (buttonstate)
		need_redraw();
}

void mousebuttonfunc(int button, int state, int x, int y)
{
	mouseClickBar = true;
	if(TwEventMouseButtonGLUT(button, state, x, y))
	{
		TwRefreshBar(bar);
		need_redraw();
		return;
	}
	mouseClickBar = false;

	static timestamp last_click_time;
	static unsigned last_click_buttonstate = 0;
	static float doubleclick_threshold = 0.4f;

	if (glutGetModifiers() & GLUT_ACTIVE_CTRL)
		buttonstate |= (1 << 30);
	else
		buttonstate &= ~(1 << 30);

	if (state == GLUT_DOWN) {
		buttonstate |= (1 << button);
		if (buttonstate == last_click_buttonstate &&
			now() - last_click_time < doubleclick_threshold) {
				doubleclick(button, x, y);
				last_click_buttonstate = 0;
		} else {
			last_click_time = now();
			last_click_buttonstate = buttonstate;
		}
	} else {
		buttonstate &= ~(1 << button);
	}
	if (button == 3 || button == 4)
		buttonstate |= (1 << button);


	mousehelperfunc(x, y);

	if (buttonstate & ((1 << 3) | (1 << 4))) // Wheel
		need_redraw();
	if ((buttonstate & 7) && (buttonstate & (1 << 30))) // Light
		need_redraw();
	if (buttonstate & ((1 << 3) | (1 << 4)))
		buttonstate &= ~((1 << 3) | (1 << 4));
}


// Idle callback
void idle()
{
	xform tmp_xf = global_xf;
	if (current_mesh >= 0)
		tmp_xf = global_xf * xforms[current_mesh];

	if (camera.autospin(tmp_xf))
		need_redraw();
	else
		usleep(10000);

	if (current_mesh >= 0) {
		xforms[current_mesh] = inv(global_xf) * tmp_xf;
		update_bsph();
	} else {
		global_xf = tmp_xf;
	}



	tmp_xf = global_xf * xforms[0];
	meshes[0]->eyePosition = inv(tmp_xf) * vec3(0.f);


}


// Keyboard
#define Ctrl (1-'a')
void keyboardfunc(unsigned char key, int, int)
{
	switch (key) {
	case ' ':
		if (current_mesh < 0)
			resetview();
		else
			set_current(-1);
		break;
		
	case 'I':
		dump_image(); break;
	case Ctrl+'x':
		save_xforms();
		break;
	case Ctrl+'c':
		save_cam_xform();
		break;
	case '\033': // Esc
		exit(0);
	case '0':
		toggle_vis(9); break;
	case '-':
		toggle_vis(10); break;
	case '=':
		toggle_vis(11); break;
	default:
		if (key >= '1' && key <= '9') {
			int m = key - '1';
			toggle_vis(m);
		}
	}
	need_redraw();
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-grab] infile...\n", myname);
	exit(1);
}

void mousePassiveMotionFunc(int x, int y)
{
	TwMouseMotion(x, y);
	TwRefreshBar(bar);
	need_redraw();

}



int _tmain(int argc, char* argv[])
{
	glutInitWindowSize(800, 800);

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInit(&argc, argv);
	




	picName.push_back("building_cross");
	picName.push_back("campus_cross");
	picName.push_back("galileo_cross");
	picName.push_back("kitchen_cross");
	picName.push_back("stpeters_cross");
	picName.push_back("uffizi_cross");
	picName.push_back("grace_cross");
	picName.push_back("rnl_cross");
	


	char *direction = {"..\\Model\\"};

	char *fileName = { "dragon" };
//	char *fileName = { "chicNN" };
//	char *fileName = { "cubeNN" };
//	char *fileName = { "headMaxNN" };
//	char *fileName = { "rocker-arm" };
//	char *fileName = { "TableclothNN" };
	

	
	

//	char *fileName = { "flubberNN" };
//	char *fileName = { "BolNHNN" };
//	char *fileName = { "cubeNN" };
//	char *fileName = { "cubeBigNN" };
//	char *fileName = { "cube rounded" };
//	char *fileName = { "TableclothNN" };
//	char *fileName = { "sphere_1W" };

//	char *fileName = { "HippoNN" };
//	char *fileName = { "maxplanckNN" };
//	char *fileName = { "screwdriverNN" };
	
	
//		char *fileName = { "Female_dancer88" };

//	char *fileName = { "794_lagomaggioreNN" };
//	char *fileName = { "circularNN" };
//	char *fileName = { "vase-lionNN" };
//	char *fileName = { "PoseidoBignNN" };
//	char *fileName = { "grog" };
	
//	char *fileName = { "794_lagomaggioreSNN" };
//	char *fileName = { "PoseidonNN" };
//	char *fileName = { "PosNN" };
//	char *fileName = { "fourNN" };
	
//	char *fileName = { "Bol1NN" };
//	char *fileName = {"ArmadilloNN"};
//	char *fileName = {"vase-lionNN"};
//	char *fileName = {"794_lagomaggioreSNN"};
//	char *fileName = { "sphereNN" };
//	char *fileName = { "sphere_1W" };
//	char *fileName = { "sphere_2W" };
//	char *fileName = { "sphere_4W" };
	
	/*

	ofstream neighbors_if("..\\Model\\BolNN.txt");

	for (int x = 0; x < 1024; x++)
	{
		for (int y = 0; y < 1024; y++)
		{

			float xx = x;
			float yy = y;

			vec2 coord(xx, yy);

			normalize(coord);
			float angle = acos(coord[1]);



			neighbors_if << angle << "	";



		}
		neighbors_if << endl;
	}

	*/




	glutCreateWindow(fileName);


	mouseXInc = 0;
	mouseYInc = 0;

	grab_only = false;

	TriMesh *themesh = TriMesh::read(direction, fileName);


	meshes.push_back(themesh);
	radianceScaling.inite(meshes[0]);
	

	/*
	Color text1(100.f, -128.f, -128.f);
	Color text2(100.f, -60.f, 128.f);
	Color text3(100.f, 128.f, -60.f);

	text1 = Color(0.f, 0.f, 0.f);
	text2 = Color(50.f, 0.f, 0.f);
	text3 = Color(100.f, 0.f, 0.f);

	text1 = text1.convert(Color::CIELAB, Color::RGB);
	text2 = text2.convert(Color::CIELAB, Color::RGB);
	text3 = text3.convert(Color::CIELAB, Color::RGB);
	
	*/

	

	if (false)
	{


		vec3 max(0.0, 0.0, 0.0);
		vec3 min(0.0, 0.0, 0.0);
		vec3 average(0.0, 0.0, 0.0);
		for (int i = 0; i < themesh->vertices.size(); i++)
		{
			average += themesh->vertices[i];
			for (int d = 0; d < 3; d++)
			{
				if (themesh->vertices[i][d] > max[d])
					max[d] = themesh->vertices[i][d];
				if (themesh->vertices[i][d] < min[d])
					min[d] = themesh->vertices[i][d];
			}
		}

		average /= themesh->vertices.size();


		min -= average;
		max -= average;
		for (int d = 0; d < 3; d++)
			cout << min[d] << " " << max[d] << " " << endl;
		/*
		themesh->vertices.clear();
		themesh->vertices.resize(radianceScaling.evnlight.directionNum);
		themesh->faces.clear();

		for (int i = 0; i < radianceScaling.evnlight.directionNum; i++)
		{
			themesh->vertices[i] = radianceScaling.evnlight.directions[i];
		}

		for (int i = 0; i < radianceScaling.evnlight.facePoints.size(); i++)
		{
			themesh->faces.push_back(trimesh::TriMesh::Face(radianceScaling.evnlight.facePoints[i][0], radianceScaling.evnlight.facePoints[i][1], radianceScaling.evnlight.facePoints[i][2]));
		}

		*/

		float a2 = 1.0;
		float b2 = 5.0;


		for (int i = 0; i < themesh->vertices.size(); i++)
		{
			vec3 temp = themesh->vertices[i];
	//		temp -= average;

			float translate = temp[1] * temp[1] * a2 + (temp[0] * temp[0] + temp[2] * temp[2])*b2;
			translate = 1.0 / sqrt(translate);
			themesh->vertices[i] = temp * translate;
		}


		/*
	
		ofstream BVH_if("..\\Model\\sdf.obj");
		BVH_if << "# Vertices: " << themesh->vertices.size() << endl;
		BVH_if << "# Faces : " << themesh->faces.size() << endl;



		for (int i = 0; i < themesh->vertices.size(); i++)
		{
			BVH_if << "v " << setiosflags(ios::fixed) << themesh->vertices[i][0] << " " << themesh->vertices[i][1] << " " << themesh->vertices[i][2] << endl;

		}


		for (int i = 0; i < themesh->faces.size(); i++)
		{
			BVH_if << "f " << themesh->faces[i][0]+1 << " " << themesh->faces[i][1] +1<< " " << themesh->faces[i][2]+1 << endl;

		}

		*/

	}
	

	themesh->need_normals();
	themesh->need_tstrips();
	themesh->need_bsphere();


	xforms.push_back(xform());
	visible.push_back(true);
	filenames.push_back(fileName);





	glutDisplayFunc(redraw);
	glutMouseFunc(mousebuttonfunc);
	glutMotionFunc(mousemotionfunc);


	glutKeyboardFunc(keyboardfunc);
	glutPassiveMotionFunc(mousePassiveMotionFunc);
	glutIdleFunc(idle);


	resetview();

	const GLubyte* name = glGetString(GL_VENDOR);
	const GLubyte* identifier = glGetString(GL_RENDERER);
	const GLubyte* version = glGetString(GL_VERSION);
	const GLubyte* gluVersion = gluGetString(GLU_VERSION);
	const GLubyte* glewVersion = glewGetString(GLEW_VERSION);
	const GLubyte* glExt = glGetString(GL_EXTENSIONS);

	printf("vendor: %s\n", name);
	printf("identifer: %s\n", identifier);
	printf("opengl version: %s\n", version);
	printf("glu version: %s\n", gluVersion);
	printf("glew version: %s\n", glewVersion);

	const GLubyte *renderer = glGetString(GL_RENDERER);
	const GLubyte *vendor = glGetString(GL_VENDOR);
	const GLubyte *glslVersion =
		glGetString(GL_SHADING_LANGUAGE_VERSION);
	GLint major, minor;
	glGetIntegerv(GL_MAJOR_VERSION, &major);
	glGetIntegerv(GL_MINOR_VERSION, &minor);
	cout << "GL Vendor    :" << vendor << endl;
	cout << "GL Renderer  : " << renderer << endl;
	cout << "GL Version (string)  : " << version << endl;
	cout << "GL Version (integer) : " << major << "." << minor << endl;
	cout << "GLSL Version : " << glslVersion << endl;

	initTwBar();

	glutMainLoop();
}




