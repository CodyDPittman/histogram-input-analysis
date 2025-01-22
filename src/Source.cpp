#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cfloat>
#include <cmath>
#include <gl/glew.h>
#include <GL/glut.h>
#include "utility.h"

/// \file
/// This is the main file for the project Input Analysis.
/// 
/// \author Cody Pittman
/// \version 1.0
/// \date 10/16/23
/// 


using namespace std;

/// \param File names
char* fileName1, * fileName2;
string currentFile, distributionType;

/// \param Input data
float* dataset; ///< Array used for X data.
int numDataPoints;
float minimum, maximum; ///< Parameters for max and min X data.

/// \param Histogram
int numIntervals = 30;
float* endPoints; ///< Array used for end points.
float* prob; ///< Array used for probability.
float maxProb = 1;

/// \param Theoretical distributions
int curveType = 0;
int numCurvePoints = 100;
float* curveX = new float[numCurvePoints]; ///< Array used for x values for PDF function.
float* curveY = new float[numCurvePoints]; ///< Array used for y values for PDF function.

/// \param Parameters
float mu = 0, sigma = 1; // Normal distribution
float lamda = 1; // Exponential distribution
float parameterStep = 0.05; // Step size for changing parameter values
float beta = 1.25;

/// Drawing parameters
int width = 800, height = 600;
float world_x_min, world_x_max, world_y_min, world_y_max;
float axis_x_min , axis_x_max, axis_y_min, axis_y_max;

///Menu variables
GLuint window, screen, command;
int fileMenu, mainMenu, histogramMenu, parameterMenu, distributionMenu;

/// \brief Compute all the points for normal distribution
void computeNormalFunc(float mu, float sigma)
{
	const double pi = 3.14159265358979323846;
	const double e = 2.71828182845904523536;
	float stepsize;

	// Delete previously allocated memory
	if (curveX != NULL)
	{
		delete[] curveX;
	}
	if (curveY != NULL)
	{
		delete[] curveY;
	}

	curveX = new float[numCurvePoints];
	curveY = new float[numCurvePoints];

	// Determine the step size and compute the arrays curveX and curveY
	// (numCurvePoints).
	stepsize = (maximum - minimum) / (numCurvePoints - 1);
	

	for (int i = 0; i < numCurvePoints; i++)
	{
		curveX[i] = minimum + (i * stepsize);
		curveY[i] = (1 / (sigma * sqrt(2 * pi))) * pow(e, -1*(pow(curveX[i] - mu, 2.0) / 2 * pow(sigma, 2.0)));
	}

}

/// Compute all the points for exponential distribution
void computeExponentialFunc(float lamda)
{
	const double e = 2.71828182845904523536;
	float stepsize;

	// Delete previously allocated memory
	if (curveX != NULL)
	{
		delete[] curveX;
	}
	if (curveY != NULL)
	{
		delete[] curveY;
	}

	curveX = new float[numCurvePoints];
	curveY = new float[numCurvePoints];


	// Determine the step size and compute the arrays curveX and curveY
	// (numCurvePoints).
	stepsize = (maximum - minimum) / (numCurvePoints - 1);

	for (int i = 0; i < numCurvePoints; i++)
	{
		curveX[i] = minimum + (i * stepsize);
		curveY[i] = (1 / beta) * pow(e, -1 * curveX[i] / beta);
	}

}

/// Main display function
void display(void)
{
	/* clear all pixels */
	glClear(GL_COLOR_BUFFER_BIT);

	// Reset Modelview matrix.
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glLineWidth(1);
	glColor3f(1, 1, 1);

	// Draw x and y axes
	glBegin(GL_LINES);
	glVertex2f(axis_x_min, 0);
	glVertex2f(axis_x_min, axis_y_max + .05); //Y-Axis 
	glVertex2f(axis_x_min, 0);
	glVertex2f(axis_x_max, 0); //X-Axis
	//glVertex2()
	glEnd();

	// Display the maximum probability value
	glRasterPos2f(axis_x_min, axis_y_max);
	printString(to_string(maxProb));


	// Draw probability histogram
	glColor3f(0.0, 1.0, 0.0);
	for (int i = 0; i < numIntervals; i++)
	{
		//glRectf(endPoints[i], prob[i], endPoints[i + 1], prob[i]);

		glBegin(GL_LINE_LOOP);
		glVertex2f(endPoints[i], 0);
		glVertex2f(endPoints[i + 1], 0);
		glVertex2f(endPoints[i + 1], prob[i]);
		glVertex2f(endPoints[i], prob[i]);
		glEnd();
	}

	// Draw the theoretical distribution using thicker lines.
	glLineWidth(2);
	glColor3f(1.0, 0.0, 0.0);
	if (curveType == 1 || curveType == 2)
	{
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < numCurvePoints; i++)
		{
			glVertex2f(curveX[i], curveY[i]);
		}
		glEnd();

	}
	
	// Compute the top-left position of the annotation
	glColor3f(1.0, 1.0, 1.0);
	glRasterPos2f(axis_x_min, axis_y_max + .1);
	printString("Probability Density");
	glRasterPos2f(axis_x_max - .2 , 0.01);
	printString("Data");

	// File Name
	glColor3f(0.0, 1.0, 0.0);
	glRasterPos2f(axis_x_max - 3, axis_y_max + .08);
	printString("File: ");
	printString(currentFile);

	// Minimum
	glRasterPos2f(axis_x_max - 3, axis_y_max + .06);
	printString("Min: ");
	printString(to_string(endPoints[0]));

	// Maximum
	glRasterPos2f(axis_x_max - 3, axis_y_max + .04);
	printString("Max: ");
	printString(to_string(endPoints[numIntervals]));

	// Number of intervals
	glRasterPos2f(axis_x_max - 3, axis_y_max + .02);
	printString("Number of Intervals: ");
	printString(to_string(numIntervals));
	
	if (currentFile == "4.dat")
	{
		glColor3f(1, 1, 1);
		glRasterPos2f(axis_x_max, 0.02);
		printString("Data");

		// File Name
		glColor3f(0.0, 1.0, 0.0);
		glRasterPos2f(axis_x_max - .2 , axis_y_max);
		printString("File: ");
		printString(currentFile);

		// Minimum
		glRasterPos2f(axis_x_max -.2, axis_y_max - .1);
		printString("Min: ");
		printString(to_string(endPoints[0]));

		// Maximum
		glRasterPos2f(axis_x_max-.2, axis_y_max -.2);
		printString("Max: ");
		printString(to_string(endPoints[numIntervals]));

		// Number of intervals
		glRasterPos2f(axis_x_max-.2, axis_y_max -.3);
		printString("Number of Intervals: ");
		printString(to_string(numIntervals));
	}


	// Draw theoretical distributions
	glColor3f(1.0, 0.0, 0.0);

	if (distributionType == "Normal")
	{
		glRasterPos2f(axis_x_max - 3, axis_y_max);
		printString("Distribution: ");
		printString(distributionType);
		glRasterPos2f(axis_x_max - 3, axis_y_max -.02);
		printString("Mu: ");
		printString(to_string(mu));
		glRasterPos2f(axis_x_max - 3, axis_y_max -.04);
		printString("Sigma: ");
		printString(to_string(sigma));
		//glutPostRedisplay();
	}
	if (distributionType == "Exponential")
	{
		glRasterPos2f(axis_x_max - 3, axis_y_max);
		printString("Distribution: ");
		printString(distributionType);
		glRasterPos2f(axis_x_max - 3, axis_y_max -.02);
		printString("Beta: ");
		printString(to_string(beta));
	}

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(world_x_min, world_x_max, world_y_min, world_y_max);

	glFlush();
	glutSwapBuffers();
}

void init(void)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
}

/// Compute the probability for the histogram (vertical axis)
void computeProbability(int numIntervals)
{
	float range, binWidth, area, height;
	float occurance;

	// Range is the maximum - minimum of the data
	range = maximum - minimum;

	// Bin width for each interval
	binWidth = range / numIntervals;


	// Delete previously allocated memory
	if (endPoints != NULL)
	{
		delete[] endPoints;
	}
	if (prob != NULL)
	{
		delete[] prob;
	}

	// Determine the end points for each interval (update the array endPoints).
	endPoints = new float[numIntervals + 1]; // two endpoints per interval so we need +1 for last interval

	for (int i = 0; i < numIntervals; i++)
	{
		// First end point always the minimum of the dataset
		if (i == 0)
			endPoints[i] = minimum;

		// Subsequent end points always the previous endpoint + bin width
		endPoints[i + 1] = endPoints[i] + binWidth;

		cout << "Endpoint " << i << ": [" << endPoints[i] << "," << endPoints[i + 1] << "]" << endl;
	}

	// Compute the probability for each interval (update the array prob).
	prob = new float[numIntervals];
	
	// Determine the amount of data points that fall in each interval
	for (int i = 0; i < numIntervals; i++)
	{
		occurance = 0;
		for (int j = 0; j < numDataPoints; j++)
		{
			if (dataset[j] > endPoints[i] && dataset[j] < endPoints[i + 1])
				occurance++;
		}
		prob[i] = occurance / numDataPoints; // the area of each rectangle

		//cout << "Probability of interval " << i << " is " << prob[i] << endl;
		//cout << "A total of " << occurance << " datapoints lie inside of interval " << i << endl;
	}

	// Determine the height
	for (int i = 0; i < numIntervals; i++)
	{
		prob[i] = prob[i] / binWidth;
		//cout << "Probability Density of interval " << i << " is " << prob[i] << endl;
	}

	float greatestProbDens = 0.0;

	for (int i = 0; i < numIntervals; i++)
	{
		if (prob[i] > greatestProbDens)
			greatestProbDens = prob[i];
	}

	maxProb = greatestProbDens;
	axis_y_max = maxProb;
	world_y_max = axis_y_max + .15;
}

/// Function for reading a file.
void readFile(string fileName)
{
	ifstream inFile(fileName);
	if (!inFile.is_open())
	{
		cout << fileName << " couldn't be opened.\n";
		system("pause");
		exit(1);
	}
	inFile >> numDataPoints;

	// Memory management.
	if (dataset != NULL)
		delete dataset;
	dataset = new float[numDataPoints];
	minimum = FLT_MAX;
	maximum = -FLT_MAX;

	// Read the data and compute minimum and maximum.
	for (int i = 0; i < numDataPoints; i++)
	{
		inFile >> dataset[i];
		if (dataset[i] < minimum)
			minimum = dataset[i];
		if (dataset[i] > maximum)
			maximum = dataset[i];
		//cout << dataset[i] << endl;
	}
	
	// Compute the limits for the axes and world.
	axis_x_max = maximum + (maximum - minimum) * .005;
	axis_x_min = minimum - (maximum - minimum) * .005;

	world_x_max = maximum + (maximum - minimum) * .05;
	world_x_min = minimum - (maximum - minimum) * .05;


	axis_y_min = 0; // minimum probability is always 0 (0%)
	world_y_min = axis_y_min - .05;

	// Compute the histogram
	computeProbability(numIntervals);

	currentFile = fileName;
}

void keyboard(unsigned char key, int x, int y)
{

	if (key == 'q' || key == 'Q' || key == 27)
		exit(0);
}

/// <summary>
/// Function for accessing arrow key functions.
/// </summary>
/// <param name="id"></param>
void specialKey(int key, int x, int y) // for the arrow keys
{
	// Increasing the parameters based upon the step size chosen

	switch (key)
	{
	case GLUT_KEY_RIGHT:
		//cout << "Input detected" << endl;
		if (distributionType == "Normal")
		{
			// Increase mu
			mu += parameterStep;

			if (mu > 5)
				mu = 5;
		}
		break;
	case GLUT_KEY_LEFT:
		//cout << "Input detected" << endl;
		if (distributionType == "Normal")
		{
			// Decrease mu
			mu -= parameterStep;

			if (mu < 0)
				mu = 0;
		}
		break;
	case GLUT_KEY_UP:
		//cout << "Input detected" << endl;
		if (distributionType == "Normal")
		{
			// Increase sigma
			sigma += parameterStep;

			if (sigma > 3)
				sigma = 3;
		}
		else if (distributionType == "Exponential")
		{
			// Increase beta
			beta += parameterStep;

			if (beta > 6)
				beta = 6;
		}
		break;
	case GLUT_KEY_DOWN:
		//cout << "Input detected" << endl;
		if (distributionType == "Normal")
		{
			// Decrease sigma
			sigma -= parameterStep;

			if (sigma < .02)
				sigma = .02;
		}
		else if (distributionType == "Exponential")
		{
			// Decrease beta
			beta -= parameterStep;

			if (beta < .1)
				beta = .1;
		}
		break;
	default:
		break;

	}

	if (curveType == 1)
	{
		computeNormalFunc(mu, sigma);
	}
	if (curveType == 2)
	{
		computeExponentialFunc(beta);
	}
	
	
	// Update the parameters and theoretical distributions
	glutPostRedisplay();
}

void topMenuFunc(int id)
{
	exit(0);
}

/// Functions associated with reading each file selection in the File sub menu.
void fileMenuFunction(int id)
{
	// Read file.
	switch (id)
	{
	case 0:
		readFile("normal.dat");
		break;
	case 1:
		readFile("expo.dat");
		break;
	case 2:
		readFile("4.dat");
		break;
	case 3:
		readFile("17.dat");
		break;
	default:
		break;
	}


	// Update projection since the data has changed.
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(world_x_min, world_x_max, world_y_min, world_y_max);


	glutPostRedisplay();
}
/// Menu for distribution selection.
void funcMenuFunction(int id)
{
	// Different theoretical distributions
	switch (id)
	{
	case 0: 
		distributionType = "Normal";
		curveType = 1;
		computeNormalFunc(mu, sigma);
		break;
	case 1:
		distributionType = "Exponential";
		curveType = 2;
		computeExponentialFunc(lamda);
		break;
	default:
		break;
	}

	// Update project since the data has changed.
	glutPostRedisplay();
}
/// Menu for interval selection.
void histogramMenuFunction(int id)
{
	// Update the number of intervals and recompute the histogram.
	// Update projection since the histogram has changed due to the change of number
	// of bars.
	switch (id)
	{
	case 0:
		numIntervals = 30;
		computeProbability(numIntervals);
		break;
	case 1:
		numIntervals = 40;
		computeProbability(numIntervals);
		break;
	case 2:
		numIntervals = 50;
		computeProbability(numIntervals);
		break;
	default:
		break;
	}

	glutPostRedisplay();
}
/// Menu for parameter step.
void parameterStepMenuFunction(int id)
{
	// Update the parameter step size.
	switch (id)
	{
	case 0:
		parameterStep = 0.01;
		break;
	case 1:
		parameterStep = 0.02;
		break;
	case 2:
		parameterStep = 0.05;
	}
	glutPostRedisplay();
}

void createMenu(int value)
{
	// Create menus
	keyboard((unsigned char)value, 0, 0);

	
}
/// <summary>
/// Function for reshaping the window when resized.
/// </summary>
void reshape(int w, int h)
{
	// Set matrix mode and projection.
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (w <= h)
		glOrtho(-2.0, 2.0, -2.0 * (GLfloat)h / (GLfloat)w,
			2.0 * (GLfloat)h / (GLfloat)w, -10.0, 10.0);
	else
		glOrtho(-2.0 * (GLfloat)w / (GLfloat)h,
			2.0 * (GLfloat)w / (GLfloat)h, -2.0, 2.0, -10.0, 10.0);
	glMatrixMode(GL_MODELVIEW);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(world_x_min, world_x_max, world_y_min, world_y_max);

	glutPostRedisplay();
}

int main(int argc, char** argv)
{
	//GLUT initalization
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(width, height);

	
	//Displaying the Window
	window = glutCreateWindow("Input Analysis");
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);

	// For having an initial graph display
	//computeProbability(numIntervals);
	readFile("normal.dat");
	

	//Callbacks
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(specialKey);
	
	//File Menu Selections
	fileMenu = glutCreateMenu(fileMenuFunction);
	glutAddMenuEntry("normal.dat", 0);
	glutAddMenuEntry("expo.dat", 1);
	glutAddMenuEntry("Data File 4", 2);
	glutAddMenuEntry("Data File 17", 3);

	//Distribution Sub Menu Selections
	distributionMenu = glutCreateMenu(funcMenuFunction);
	glutAddMenuEntry("Normal", 0);
	glutAddMenuEntry("Exponential", 1);

	//Histogram Sub Menu Selections
	histogramMenu = glutCreateMenu(histogramMenuFunction);
	glutAddMenuEntry("30", 0);
	glutAddMenuEntry("40", 1);
	glutAddMenuEntry("50", 2);

	//Parameter Step Sub Menu Selections
	parameterMenu = glutCreateMenu(parameterStepMenuFunction);
	glutAddMenuEntry("0.01", 0);
	glutAddMenuEntry("0.02", 1);
	glutAddMenuEntry("0.05", 2);

	//Main Menu Selections
	mainMenu = glutCreateMenu(createMenu);
	glutAddSubMenu("Files", fileMenu);
	glutAddSubMenu("Distribution", distributionMenu);
	glutAddSubMenu("Histogram", histogramMenu);
	glutAddSubMenu("Parameter Step:", parameterMenu);
	glutAddMenuEntry("Exit", 27);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	
	init();

	//GLUT main loop
	glutMainLoop();
}