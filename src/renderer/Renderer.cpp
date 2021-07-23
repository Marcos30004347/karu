#include "renderer/Renderer.hpp"

#include <iostream>
#include <vector>
#include "linmath.h"

static const char* vertex_shader_text =
"#version 110\n"
"uniform mat4 MVP;\n"
"attribute vec2 vPos;\n"
"void main()\n"
"{\n"
"    gl_Position = MVP * vec4(vPos, 0.0, 1.0);\n"
"}\n";
 
static const char* fragment_shader_text =
"#version 110\n"
"void main()\n"
"{\n"
"    gl_FragColor = vec4(1.0, 1, 0, 1.0);\n"
"}\n";


static void error_callback(int error, const char* description)
{
		fprintf(stderr, "Error: %s\n", description);
}

static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
		if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
				glfwSetWindowShouldClose(window, GLFW_TRUE);
}

namespace karu 
{

Renderer::Renderer(size_t width, size_t heigth)
{
	window = nullptr;

	glfwSetErrorCallback(error_callback);

	if (!glfwInit())
	{
		fprintf(stderr, "GLEW initialization error\n");
		exit(EXIT_FAILURE);
	}
	
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

	this->window = glfwCreateWindow(width, heigth, "Karu", NULL, NULL);
	
	if (!this->window)
	{
		fprintf(stderr, "Window initialization error\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	
	glfwSetKeyCallback(window, keyCallback);

	glfwMakeContextCurrent(window);
	int gladInitRes = gladLoadGL();

	if (!gladInitRes) {
		fprintf(stderr, "Unable to initialize glad\n");
		glfwDestroyWindow(window);
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
}

Renderer::~Renderer()
{
	glfwDestroyWindow(window);
	glfwTerminate();
}

void Renderer::draw2dPoints(std::vector<algebra::Vector> points, float centerX, float centerY, float normX, float normY)
{
	std::vector<float> data;

	for(int i=0; i<points.size(); i++)
	{
		data.push_back((points[i][0] - centerX)/normX);
		data.push_back((points[i][1] - centerY)/normY);
	}

	glfwSwapInterval(1);

	GLuint vertex_buffer, vertex_shader, fragment_shader, program;
	GLint mvp_location, vpos_location, vcol_location;

	glGenBuffers(1, &vertex_buffer);

	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
	glBufferData(GL_ARRAY_BUFFER, data.size()*sizeof(float), data.data(), GL_STATIC_DRAW);

	vertex_shader = glCreateShader(GL_VERTEX_SHADER);

	glShaderSource(vertex_shader, 1, &vertex_shader_text, NULL);
	glCompileShader(vertex_shader);

	fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);

	glShaderSource(fragment_shader, 1, &fragment_shader_text, NULL);
	glCompileShader(fragment_shader);

	program = glCreateProgram();

	glAttachShader(program, vertex_shader);
	glAttachShader(program, fragment_shader);
	
	glLinkProgram(program);

	mvp_location = glGetUniformLocation(program, "MVP");
	vpos_location = glGetAttribLocation(program, "vPos");

	glEnableVertexAttribArray(vpos_location);
	glVertexAttribPointer(vpos_location, 2, GL_FLOAT, GL_FALSE, sizeof(float)*2, (void*) 0);

	glClearColor(0.15f, 0.15f, 0.15f, 1.0f);


	mat4x4 m, p, mvp;
	float ratio;
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	ratio = width / (float) height;
	mat4x4_identity(m);
	mat4x4_ortho(p, -ratio, ratio, -1.f, 1.f, 1.f, -1.f);
	mat4x4_mul(mvp, p, m);

	glPointSize(2);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	while (!glfwWindowShouldClose(window))
	{

			glViewport(0, 0, width, height);
			glClear(GL_COLOR_BUFFER_BIT);
		

			glUseProgram(program);
			glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (const GLfloat*) mvp);
			
			glDrawArrays(GL_POINTS, 0, points.size());

			glfwSwapBuffers(window);
			glfwPollEvents();
	}
}

}
