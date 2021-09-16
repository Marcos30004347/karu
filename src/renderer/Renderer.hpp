#pragma once

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <GL/glu.h>
#include <vector>

#include "algebra/matrix/Matrix.hpp"
#include "camera/Camera.hpp"


namespace karu
{
	


class Renderer
{
	GLFWwindow* window;
	u64 width;
	u64 heigth;
public:
	Renderer(size_t width, size_t heigth);
	~Renderer();
	void draw2dPoints(Camera& cam, std::vector<algebra::Matrix> points);
};

} // namespace karu
