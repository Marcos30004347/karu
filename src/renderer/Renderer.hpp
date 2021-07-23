

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <GL/glu.h>
#include <vector>

#include "algebra/vector/Vector.hpp"



namespace karu
{
	


class Renderer
{
	GLFWwindow* window;

public:
	Renderer(size_t width, size_t heigth);
	~Renderer();
	void draw2dPoints(std::vector<algebra::Vector> points, float centerX = 0, float centerY = 0, float normX = 1, float normY = 1);
};

} // namespace karu
