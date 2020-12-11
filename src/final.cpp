#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include "glm/ext.hpp"
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <stdlib.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>

#define WIDTH 500
#define HEIGHT 500

std::array<std::array<float, HEIGHT>, WIDTH> depthBuffer;

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	std::vector<float> v;
	if (numberOfValues > 1) {
		float d = (to - from) / (numberOfValues - 1);
		for (size_t i = 0; i < numberOfValues; i++) v.push_back(from + (i * d));
	} else {
		v.push_back(from);
	}
	return v;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<float> x;
	std::vector<float> y;
	std::vector<float> z;

	x = interpolateSingleFloats(from.x, to.x, numberOfValues);
	y = interpolateSingleFloats(from.y, to.y, numberOfValues);
	z = interpolateSingleFloats(from.z, to.z, numberOfValues);

	std::vector<glm::vec3> v;
	for (size_t i = 0; i < numberOfValues; i++) {
		v.push_back(glm::vec3(x.at(i), y.at(i), z.at(i)));
	}

	return v;
}

bool inRange(float n, float x0, float x1) {
	return (n >= fmin(x0, x1) && n < fmax(x0, x1));
}

void drawLine(DrawingWindow &window, CanvasPoint::CanvasPoint to, CanvasPoint::CanvasPoint from, Colour::Colour colour) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberofSteps = fmax(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/numberofSteps;
	float yStepSize = yDiff/numberofSteps;

	uint32_t intColour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
	std::vector<float> ds = interpolateSingleFloats(from.depth, to.depth, numberofSteps);
	for (size_t i = 0; i < numberofSteps; i++) {
		float x = from.x + (xStepSize * i);
		float y = from.y + (yStepSize * i);
		float d = -1/ds.at(i);

		if (inRange(x, 0, WIDTH-1) && inRange(y, 0, HEIGHT-1) && d > depthBuffer[x][y]) {
			depthBuffer[x][y] = d;
			window.setPixelColour(round(x), round(y), intColour);
		}
	}
}

void drawUnfilledTriangle(DrawingWindow &window, CanvasTriangle::CanvasTriangle triangle, Colour::Colour colour) {
	drawLine(window, triangle.vertices[0], triangle.vertices[1], colour);
	drawLine(window, triangle.vertices[1], triangle.vertices[2], colour);
	drawLine(window, triangle.vertices[2], triangle.vertices[0], colour);
}

std::array<CanvasPoint, 3> orderTriangleVertices(CanvasTriangle::CanvasTriangle triangle) {
	if (triangle.vertices[0].y < triangle.vertices[1].y) {
		if (triangle.vertices[2].y < triangle.vertices[0].y) std::swap(triangle.vertices[0], triangle.vertices[2]);
	} else if (triangle.vertices[1].y < triangle.vertices[2].y) {
		std::swap(triangle.vertices[0], triangle.vertices[1]);
	} else std::swap(triangle.vertices[0], triangle.vertices[2]);
	if (triangle.vertices[2].y < triangle.vertices[1].y) std::swap(triangle.vertices[1], triangle.vertices[2]);
	return triangle.vertices;
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle::CanvasTriangle triangle, Colour::Colour colour) {
	triangle.vertices = orderTriangleVertices(triangle);

	// Calculate x coordinate of 4th point
	float m = (triangle.vertices[2].x - triangle.vertices[0].x) / (triangle.vertices[2].y - triangle.vertices[0].y);
	float x = triangle.vertices[0].x + (m * (triangle.vertices[1].y - triangle.vertices[0].y));

	// Calculate depth of 4th point
	m = (triangle.vertices[2].depth - triangle.vertices[0].depth) / (triangle.vertices[2].y - triangle.vertices[0].y);
	float depth = triangle.vertices[0].depth + (m * (triangle.vertices[1].y - triangle.vertices[0].y));

	CanvasPoint::CanvasPoint point = CanvasPoint::CanvasPoint(x, triangle.vertices[1].y, depth);

	float topDY = (round(triangle.vertices[1].y) - round(triangle.vertices[0].y))+1;
	float bottomDY = (round(triangle.vertices[2].y) - round(triangle.vertices[1].y))+1;

	std::vector<float> leftLine;
	std::vector<float> leftLineDepth;
	std::vector<float> rightLine;
	std::vector<float> rightLineDepth;
	CanvasPoint::CanvasPoint leftPoint;
	CanvasPoint::CanvasPoint rightPoint;

	// Draw top half of triangle
	leftLine = interpolateSingleFloats(triangle.vertices[0].x, point.x, topDY);
	leftLineDepth = interpolateSingleFloats(triangle.vertices[0].depth, point.depth, topDY);
	rightLine = interpolateSingleFloats(triangle.vertices[0].x, triangle.vertices[1].x, topDY);
	rightLineDepth = interpolateSingleFloats(triangle.vertices[0].depth, triangle.vertices[1].depth, topDY);
	for (size_t i = 0; i < topDY; i++) {
		if (leftLine[i] > rightLine[i]) {
			leftLine[i] = ceil(leftLine[i])+1;
			rightLine[i] = floor(rightLine[i])-1;
		} else {
			leftLine[i] = floor(leftLine[i])-1;
			rightLine[i] = ceil(rightLine[i])+1;
		}
		leftPoint = CanvasPoint::CanvasPoint(leftLine[i], floor(triangle.vertices[0].y)+i, leftLineDepth[i]);
		rightPoint = CanvasPoint::CanvasPoint(rightLine[i], floor(triangle.vertices[0].y)+i, rightLineDepth[i]);
		drawLine(window, leftPoint, rightPoint, colour);
	}

	// Draw bottom half of triangle
	leftLine = interpolateSingleFloats(point.x, triangle.vertices[2].x, bottomDY);
	leftLineDepth = interpolateSingleFloats(point.depth, triangle.vertices[2].depth, bottomDY);
	rightLine = interpolateSingleFloats(triangle.vertices[1].x, triangle.vertices[2].x, bottomDY);
	rightLineDepth = interpolateSingleFloats(triangle.vertices[1].depth, triangle.vertices[2].depth, bottomDY);
	for (size_t i = 0; i < bottomDY; i++) {
		if (leftLine[i] > rightLine[i]) {
			leftLine[i] = ceil(leftLine[i])+1;
			rightLine[i] = floor(rightLine[i])-1;
		} else {
			leftLine[i] = floor(leftLine[i])-1;
			rightLine[i] = ceil(rightLine[i])+1;
		}
		leftPoint = CanvasPoint::CanvasPoint(leftLine[i], floor(point.y)+i, leftLineDepth[i]);
		rightPoint = CanvasPoint::CanvasPoint(rightLine[i], floor(triangle.vertices[1].y)+i, rightLineDepth[i]);
		drawLine(window, leftPoint, rightPoint, colour);
	}

	//drawUnfilledTriangle(window, triangle, Colour::Colour(255, 255, 255));
}

int indexTexture(TextureMap::TextureMap texture, int x, int y) {
	return x + (texture.width*y);
}

void drawTexturedTriangle(DrawingWindow &window, CanvasTriangle::CanvasTriangle triangle, TextureMap::TextureMap texture) {
	triangle.vertices = orderTriangleVertices(triangle);

	// Calculate x coordinate of 4th point
	float m = (triangle.vertices[2].x - triangle.vertices[0].x) / (triangle.vertices[2].y - triangle.vertices[0].y);
	float x = triangle.vertices[0].x + (m * (triangle.vertices[1].y - triangle.vertices[0].y));

	// Calculate depth of 4th point
	m = (triangle.vertices[2].depth - triangle.vertices[0].depth) / (triangle.vertices[2].y - triangle.vertices[0].y);
	float depth = triangle.vertices[0].depth + (m * (triangle.vertices[1].y - triangle.vertices[0].y));

	CanvasPoint::CanvasPoint point = CanvasPoint::CanvasPoint(x, triangle.vertices[1].y, depth);

	float xRatio = (point.x - triangle.vertices[0].x) / (triangle.vertices[2].x - triangle.vertices[0].x);
	float yRatio = (point.y - triangle.vertices[0].y) / (triangle.vertices[2].y - triangle.vertices[0].y);

	// Location of point on texture map
	point.texturePoint.x = triangle.vertices[0].texturePoint.x + (xRatio * (triangle.vertices[2].texturePoint.x - triangle.vertices[0].texturePoint.x));
	point.texturePoint.y = triangle.vertices[0].texturePoint.y + (yRatio * (triangle.vertices[2].texturePoint.y - triangle.vertices[0].texturePoint.y));

	float topDY    = round(triangle.vertices[1].y) - round(triangle.vertices[0].y)+1;
	float bottomDY = round(triangle.vertices[2].y) - round(triangle.vertices[1].y)+1;

	std::vector<glm::vec3> leftLine;
	std::vector<glm::vec3> rightLine;
	CanvasPoint::CanvasPoint leftPoint;
	CanvasPoint::CanvasPoint rightPoint;
	std::vector<float> xPoints;
	std::vector<float> yPoints;

	std::vector<float> leftLineDepth;
	std::vector<float> rightLineDepth;

	glm::vec3 v1(triangle.vertices[0].x, triangle.vertices[0].texturePoint.x, triangle.vertices[0].texturePoint.y);
	glm::vec3 v2(triangle.vertices[1].x, triangle.vertices[1].texturePoint.x, triangle.vertices[1].texturePoint.y);
	glm::vec3 v3(triangle.vertices[2].x, triangle.vertices[2].texturePoint.x, triangle.vertices[2].texturePoint.y);
	glm::vec3 vp(point.x, point.texturePoint.x, point.texturePoint.y);

	// Draw top half of triangle
	leftLine  = interpolateThreeElementValues(v1, v2, topDY);
	leftLineDepth = interpolateSingleFloats(triangle.vertices[0].depth, triangle.vertices[1].depth, topDY);
	rightLine = interpolateThreeElementValues(v1, vp, topDY);
	rightLineDepth = interpolateSingleFloats(triangle.vertices[0].depth, point.depth, topDY);
	for (size_t i = 0; i < topDY; i++) {
		int y = round(triangle.vertices[0].y)+i;

		leftPoint = CanvasPoint::CanvasPoint(round(leftLine[i].x), y, leftLineDepth[i]);
		leftPoint.texturePoint.x = round(leftLine[i].y);
		leftPoint.texturePoint.y = round(leftLine[i].z);

		rightPoint = CanvasPoint::CanvasPoint(round(rightLine[i].x), y, rightLineDepth[i]);
		rightPoint.texturePoint.x = round(rightLine[i].y);
		rightPoint.texturePoint.y = round(rightLine[i].z);

		int xDiff = abs(round(rightPoint.x) - round(leftPoint.x));

		std::vector<float> ds;
		if (leftPoint.x < rightPoint.x) {
			xPoints = interpolateSingleFloats(leftPoint.texturePoint.x, rightPoint.texturePoint.x, xDiff);
			yPoints = interpolateSingleFloats(leftPoint.texturePoint.y, rightPoint.texturePoint.y, xDiff);
			ds = interpolateSingleFloats(leftPoint.depth, rightPoint.depth, xDiff);
		} else {
			xPoints = interpolateSingleFloats(rightPoint.texturePoint.x, leftPoint.texturePoint.x, xDiff);
			yPoints = interpolateSingleFloats(rightPoint.texturePoint.y, leftPoint.texturePoint.y, xDiff);
			ds = interpolateSingleFloats(rightPoint.depth, leftPoint.depth, xDiff);
		}

		int canvasX;
		float d;
		for (size_t j = 0; j < xDiff; j++) {
			canvasX = round(fmin(leftPoint.x, rightPoint.x))+j;
			d = -1/ds.at(j);
			if (inRange(canvasX, 0, WIDTH-1) && inRange(y, 0, HEIGHT-1) && d > depthBuffer[canvasX][y]) {
				depthBuffer[canvasX][y] = d;
				window.setPixelColour(canvasX, y, texture.pixels[indexTexture(texture, xPoints[j], yPoints[j])]);
			}
		}
	}

	// Draw bottom half of triangle
	leftLine  = interpolateThreeElementValues(v2, v3, bottomDY);
	leftLineDepth = interpolateSingleFloats(triangle.vertices[1].depth, triangle.vertices[2].depth, bottomDY);
	rightLine = interpolateThreeElementValues(vp, v3, bottomDY);
	rightLineDepth = interpolateSingleFloats(point.depth, triangle.vertices[2].depth, bottomDY);
	for (size_t i = 0; i < bottomDY; i++) {
		int y = round(triangle.vertices[1].y)+i;

		leftPoint = CanvasPoint::CanvasPoint(round(leftLine[i].x), y, leftLineDepth[i]);
		leftPoint.texturePoint.x = round(leftLine[i].y);
		leftPoint.texturePoint.y = round(leftLine[i].z);

		rightPoint = CanvasPoint::CanvasPoint(round(rightLine[i].x), y, rightLineDepth[i]);
		rightPoint.texturePoint.x = round(rightLine[i].y);
		rightPoint.texturePoint.y = round(rightLine[i].z);

		int xDiff = abs(round(rightPoint.x) - round(leftPoint.x));

		std::vector<float> ds;
		if (leftPoint.x < rightPoint.x) {
			xPoints = interpolateSingleFloats(leftPoint.texturePoint.x, rightPoint.texturePoint.x, xDiff);
			yPoints = interpolateSingleFloats(leftPoint.texturePoint.y, rightPoint.texturePoint.y, xDiff);
			ds = interpolateSingleFloats(leftPoint.depth, rightPoint.depth, xDiff);
		} else {
			xPoints = interpolateSingleFloats(rightPoint.texturePoint.x, leftPoint.texturePoint.x, xDiff);
			yPoints = interpolateSingleFloats(rightPoint.texturePoint.y, leftPoint.texturePoint.y, xDiff);
			ds = interpolateSingleFloats(rightPoint.depth, leftPoint.depth, xDiff);
		}

		int canvasX;
		float d;
		for (size_t j = 0; j < xDiff; j++) {
			canvasX = round(fmin(leftPoint.x, rightPoint.x))+j;
			d = -1/ds.at(j);
			if (inRange(canvasX, 0, WIDTH-1) && inRange(y, 0, HEIGHT-1) && d > depthBuffer[canvasX][y]) {
				depthBuffer[canvasX][y] = d;
				window.setPixelColour(canvasX, y, texture.pixels[indexTexture(texture, xPoints[j], yPoints[j])]);
			}
		}
	}

	//drawUnfilledTriangle(window, triangle, Colour::Colour(255, 255, 255));
}

Colour colourFromName(const std::string &name, std::vector<Colour::Colour> colours) {
	for (int i = 0; i < colours.size(); i++) {
		std::vector<std::string> colourData = split(colours.at(i).name, ':');
		std::string              colourName = colourData.at(0);
		if (colourName.compare(name) == 0) {
			return colours.at(i);
		}
	}

	return Colour::Colour("Black", 0, 0, 0);
}

std::vector<ModelTriangle::ModelTriangle> readOBJ(const std::string &filename, float scale, std::vector<Colour::Colour> colours) {
	std::ifstream inputStream(filename, std::ifstream::in);
	std::string nextLine;
	std::getline(inputStream, nextLine);

	std::vector<std::string> currentLine;
	std::vector<glm::vec3> v;
	std::vector<TexturePoint::TexturePoint> vt;
	std::vector<ModelTriangle::ModelTriangle> triangles;
	Colour::Colour colour;
	std::vector<std::string> colourData;
	TextureMap::TextureMap texture;
	bool textured = false;
	std::array<std::vector<std::string>, 3>  vertices;

	bool endOfFile = false;
	while (!endOfFile) {
		endOfFile = inputStream.peek() == EOF;
		if (nextLine != "") {
			currentLine = split(nextLine, ' ');
			if (currentLine.at(0).compare("usemtl") == 0) {
				colour = colourFromName(currentLine.at(1), colours);
				if (colour.red == -1) {
					texture = TextureMap::TextureMap("models/"+split(colour.name, ':').at(1));
					textured = true;
			 	} else textured = false;
			} else if (currentLine.at(0).compare("v") == 0) {
				v.push_back(glm::vec3(std::stof(currentLine.at(1))*scale, std::stof(currentLine.at(2))*scale, std::stof(currentLine.at(3))*scale));
			} else if (currentLine.at(0).compare("vt") == 0) {
				vt.push_back(TexturePoint::TexturePoint(std::stof(currentLine.at(1))*texture.width, std::stof(currentLine.at(2))*texture.height));
			} else if (currentLine.at(0).compare("f") == 0) {
				for (int i = 0; i < 3; i++) {
					vertices[i] = split(currentLine.at(i+1), '/');
				}
				glm::vec3 v1 = v.at(std::stoi(vertices[0].at(0))-1);
				glm::vec3 v2 = v.at(std::stoi(vertices[1].at(0))-1);
				glm::vec3 v3 = v.at(std::stoi(vertices[2].at(0))-1);

				triangles.push_back(ModelTriangle::ModelTriangle(v1, v2, v3, colour));
				triangles.back().normal = glm::cross(v2-v1, v3-v1);

				if (textured) {
					triangles.back().texturePoints[0] = vt.at(std::stoi(vertices[0].at(1))-1);
					triangles.back().texturePoints[1] = vt.at(std::stoi(vertices[1].at(1))-1);
					triangles.back().texturePoints[2] = vt.at(std::stoi(vertices[2].at(1))-1);
				}
			}
		}
		std::getline(inputStream, nextLine);
	}

	inputStream.close();

	return triangles;
}

std::vector<Colour::Colour> readMTL(const std::string &filename) {
	std::ifstream inputStream(filename, std::ifstream::in);
	std::string nextLine;
	std::getline(inputStream, nextLine);

	std::vector<std::string> currentLine;
	std::vector<Colour::Colour> colours;

	bool endOfFile = false;
	while (!endOfFile) {
		endOfFile = inputStream.peek() == EOF;
		if (nextLine != "") {
			currentLine = split(nextLine, ' ');
			if (currentLine.at(0).compare("newmtl") == 0) {
				std::string name = currentLine.at(1);
				std::getline(inputStream, nextLine);
				currentLine = split(nextLine, ' ');

				// Check if texture
				if (inputStream.peek() != EOF) {
					std::getline(inputStream, nextLine);
					if (nextLine == "") {
						colours.push_back(Colour::Colour(name, std::stof(currentLine.at(1))*255, std::stof(currentLine.at(2))*255, std::stof(currentLine.at(3))*255));
					} else {
						currentLine = split(nextLine, ' ');
						// Return texture as colour with name "NAME:FILENAME" and rgb values all -1;
						colours.push_back(Colour::Colour(name+":"+currentLine.at(1), -1, -1, -1));
					}
				} else {
					colours.push_back(Colour::Colour(name, std::stof(currentLine.at(1))*255, std::stof(currentLine.at(2))*255, std::stof(currentLine.at(3))*255));
				}
			}
		} else {
			std::getline(inputStream, nextLine);
		}
	}

	inputStream.close();

	return colours;
}

glm::vec3 projectVertex(glm::vec3 vertex, glm::vec3 cameraPos, glm::mat3 cameraAngle, float f, float scale) {

	glm::vec3 camToVertex = vertex - cameraPos;

	glm::vec3 adjustedVector = camToVertex*cameraAngle;

	glm::vec3 newPos;
	newPos.x = (scale * f * adjustedVector.x / adjustedVector.z) + (WIDTH / 2);
	newPos.y = (scale * f * adjustedVector.y / adjustedVector.z) + (HEIGHT / 2);
	newPos.z = adjustedVector.z;

	return newPos;
}

CanvasTriangle projectTriangle(ModelTriangle::ModelTriangle triangle, glm::vec3 cameraPos, glm::mat3 cameraAngle, float f, float scale) {
	glm::vec3 v1 = projectVertex(triangle.vertices[0], cameraPos, cameraAngle, f, scale);
	glm::vec3 v2 = projectVertex(triangle.vertices[1], cameraPos, cameraAngle, f, scale);
	glm::vec3 v3 = projectVertex(triangle.vertices[2], cameraPos, cameraAngle, f, scale);

	CanvasPoint::CanvasPoint cp1 = CanvasPoint::CanvasPoint(v1.x, v1.y, v1.z);
	CanvasPoint::CanvasPoint cp2 = CanvasPoint::CanvasPoint(v2.x, v2.y, v2.z);
	CanvasPoint::CanvasPoint cp3 = CanvasPoint::CanvasPoint(v3.x, v3.y, v3.z);

	if (triangle.colour.red == -1) {
		cp1.texturePoint = TexturePoint::TexturePoint(triangle.texturePoints[0].x, triangle.texturePoints[0].y);
		cp2.texturePoint = TexturePoint::TexturePoint(triangle.texturePoints[1].x, triangle.texturePoints[1].y);
		cp3.texturePoint = TexturePoint::TexturePoint(triangle.texturePoints[2].x, triangle.texturePoints[2].y);
	}

	return CanvasTriangle::CanvasTriangle(cp1, cp2, cp3);
}

void drawTriangle(DrawingWindow &window, ModelTriangle::ModelTriangle triangle, CanvasTriangle::CanvasTriangle t) {
	if (triangle.colour.red == -1) {
		//drawFilledTriangle(window, t, Colour::Colour(0, 255, 0));
		TextureMap::TextureMap texture = TextureMap::TextureMap("models/"+split(triangle.colour.name, ':').at(1));
		drawTexturedTriangle(window, t, texture);
	} else drawFilledTriangle(window, t, triangle.colour);
}

RayTriangleIntersection getClosestIntersection(glm::vec3 startPos, glm::vec3 rayDirection, std::vector<ModelTriangle::ModelTriangle> triangles, float distance, int triangleIndex) {
	glm::vec3 solution;
	int solutionIndex = -1;
	for (int i = 0; i < triangles.size(); i++) {
		glm::vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
		glm::vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
		glm::vec3 SPVector = startPos - triangles[i].vertices[0];
		glm::mat3 DEMatrix(-glm::normalize(rayDirection), e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

		if ((inRange(possibleSolution.y, 0, 1) && inRange(possibleSolution.z, 0, 1) && possibleSolution.y+possibleSolution.z <= 1)
		 && (solutionIndex == -1 || possibleSolution.x < solution.x) && (distance == -1 || inRange(possibleSolution.x, 0, distance)) && i != triangleIndex) {
			solution = glm::vec3(possibleSolution.x, possibleSolution.y, possibleSolution.z);
			solutionIndex = i;
		}
	}

	if (solutionIndex != -1) {
		glm::vec3 e0 = triangles[solutionIndex].vertices[1] - triangles[solutionIndex].vertices[0];
		glm::vec3 e1 = triangles[solutionIndex].vertices[2] - triangles[solutionIndex].vertices[0];
		glm::vec3 point = triangles[solutionIndex].vertices[0] + (solution.y*e0) + (solution.z*e1);

		return RayTriangleIntersection::RayTriangleIntersection(point, solution.x, triangles[solutionIndex], solutionIndex);
	} else {
		glm::vec3 v = glm::vec3(-1,-1,-1);
		return RayTriangleIntersection::RayTriangleIntersection(v, -1, ModelTriangle::ModelTriangle(v,v,v,Colour::Colour(0,0,0)), -1);
	}
}

void drawRayTraced(DrawingWindow &window, glm::vec3 cameraPos, glm::mat3 cameraAngle, float focalLength, float scale, std::vector<ModelTriangle::ModelTriangle> triangles, glm::vec3 light) {
	window.clearPixels();

	for (int x = 0; x < WIDTH; x++) {
		for (int y = 0; y < HEIGHT; y++) {
			glm::vec3 right   = (x -  (WIDTH/2)) * cameraAngle[0] / scale;
			glm::vec3 up      = (y - (HEIGHT/2)) * cameraAngle[1] / scale;
			glm::vec3 forward = cameraAngle[2];

			glm::vec3 ray = -(right+up+forward);

			RayTriangleIntersection::RayTriangleIntersection result = getClosestIntersection(cameraPos, ray, triangles, -1, -1);

			if (result.triangleIndex != -1) {
				glm::vec3 shadowRay = light - result.intersectionPoint;
				float rayLength = glm::length(shadowRay);

				float brightness = 1/(4*M_PI*rayLength*rayLength);

				float aoi = glm::dot(glm::normalize(shadowRay), glm::normalize(result.intersectedTriangle.normal));
				aoi = 1+(aoi*0.5);
				brightness *= aoi;

				RayTriangleIntersection::RayTriangleIntersection shadowIntersect = getClosestIntersection(result.intersectionPoint, shadowRay, triangles, rayLength, result.triangleIndex);
				if (shadowIntersect.triangleIndex != -1) {
					brightness = 0;
				}

				brightness += 0.2; // Universal supplement
				if (brightness > 1) {
					brightness = 1;
				} else if (brightness < 0) {
					brightness = 0;
				}

				/*
				glm::vec3 r = glm::normalize(shadowRay) - (2*glm::normalize(result.intersectedTriangle.normal)*glm::dot(glm::normalize(shadowRay), glm::normalize(result.intersectedTriangle.normal)));
				float spec = glm::dot(ray, r);
				if (spec < 0) spec = 0;
				spec = pow(spec, 256)*255;
				if (spec > 255) spec = 255;*/
				float spec = 0;

				Colour::Colour colour = result.intersectedTriangle.colour;
				if (result.intersectedTriangle.colour.red == -1) {
					colour = Colour::Colour(0, 255, 0);
				}

				//glm::vec3 rgb = glm::vec3((255*spec)+(colour.red*(1-spec)), (255*spec)+(colour.green*(1-spec)), (255*spec)+(colour.blue*(1-spec)));
				uint32_t intColour = (255 << 24) + ((int)(fmax(colour.red,spec)*brightness) << 16) + ((int)(fmax(colour.green,spec)*brightness) << 8) + (int)(fmax(colour.blue,spec)*brightness);
				window.setPixelColour(x, y, intColour);
			}
		}
	}
}

void drawRasterized(DrawingWindow &window, glm::vec3 cameraPos, glm::mat3 cameraAngle, float focalLength, float scale, std::vector<ModelTriangle::ModelTriangle> triangles) {
	window.clearPixels();
	for (uint32_t x = 0; x < WIDTH; x++) {
		for (uint32_t y = 0; y < HEIGHT; y++) {
			depthBuffer[x][y] = 0;
		}
	}

	CanvasTriangle::CanvasTriangle t;
	for (int i = 0; i < triangles.size(); i++) {

		t = projectTriangle(triangles[i], cameraPos, cameraAngle, focalLength, scale);
		drawTriangle(window, triangles[i], t);
	}
}

void drawStroked(DrawingWindow &window, glm::vec3 cameraPos, glm::mat3 cameraAngle, float focalLength, float scale, std::vector<ModelTriangle::ModelTriangle> triangles) {
	window.clearPixels();
	for (uint32_t x = 0; x < WIDTH; x++) {
		for (uint32_t y = 0; y < HEIGHT; y++) {
			depthBuffer[x][y] = 0;
		}
	}

	CanvasTriangle::CanvasTriangle t;
	for (int i = 0; i < triangles.size(); i++) {
		t = projectTriangle(triangles[i], cameraPos, cameraAngle, focalLength, scale);
		drawUnfilledTriangle(window, t, Colour::Colour(255, 255, 255));
	}
}

void draw(DrawingWindow &window, int mode, glm::vec3 cameraPos, glm::mat3 cameraAngle, float focalLength, float scale, std::vector<ModelTriangle::ModelTriangle> triangles, glm::vec3 light) {
	if (mode == 0) drawStroked(window, cameraPos, cameraAngle, focalLength, scale, triangles);
	else if (mode == 1) drawRasterized(window, cameraPos, cameraAngle, focalLength, scale, triangles);
	else if (mode == 2) drawRayTraced(window, cameraPos, cameraAngle, focalLength, scale, triangles, light);
}

glm::mat3 lookAt(glm::vec3 cameraPos, glm::vec3 centre) {
	glm::vec3 abitrary = glm::vec3(0, -1, 0);

	glm::vec3 forward = glm::normalize(cameraPos - centre);
	glm::vec3 right = glm::normalize(glm::cross(abitrary, forward));
	glm::vec3 up = glm::normalize(glm::cross(right, forward));

	return glm::mat3(right, up, forward);
}

bool gotoPoint(glm::vec3 &cameraPos, glm::vec3 point, float speed) {
	if (glm::length(point - cameraPos) > speed) {
		cameraPos = cameraPos + (glm::normalize(point - cameraPos) * speed);
		return false;
	} else {
		cameraPos = point;
		return true;
	}
}

void update(DrawingWindow &window, glm::vec3 &cameraPos, glm::mat3 &cameraAngle, glm::vec3 &translation, glm::vec3 &rotation, glm::vec3 &orientation, std::vector<glm::vec3> cameraPath, int &pathIndex) {
	//if (gotoPoint(cameraPos, cameraPath.at(pathIndex), 0.1)) pathIndex = (pathIndex+1)%cameraPath.size();

	glm::mat3 rotMatrix;

	cameraPos = cameraPos + translation;

	if (rotation.x != 0) {
		rotMatrix = glm::mat3(
			1,                0,               0,
			0,  cos(rotation.x), sin(rotation.x),
			0, -sin(rotation.x), cos(rotation.x)
		);
		cameraPos = rotMatrix * cameraPos;
	}

	if (rotation.y != 0) {
		rotMatrix = glm::mat3(
			cos(rotation.y), 0, -sin(rotation.y),
					 0,      1,                0,
			sin(rotation.y), 0,  cos(rotation.y)
		);
		cameraPos = rotMatrix * cameraPos;
	}

	if (rotation.z != 0) {
		rotMatrix = glm::mat3(
			 cos(rotation.z), sin(rotation.z), 0,
			-sin(rotation.z), cos(rotation.z), 0,
			               0,               0, 1
		);
		cameraPos = rotMatrix * cameraPos;
	}

	if (orientation.x != 0) {
		rotMatrix = glm::mat3(
			1,                   0,                  0,
			0,  cos(orientation.x), sin(orientation.x),
			0, -sin(orientation.x), cos(orientation.x)
		);
		cameraAngle = cameraAngle * rotMatrix;
	}

	if (orientation.y != 0) {
		rotMatrix = glm::mat3(
			cos(orientation.y), 0, -sin(orientation.y),
					         0, 1,                   0,
			sin(orientation.y), 0,  cos(orientation.y)
		);
		cameraAngle = cameraAngle * rotMatrix;
	}

	if (orientation.z != 0) {
		rotMatrix = glm::mat3(
			 cos(orientation.z), sin(orientation.z), 0,
			-sin(orientation.z), cos(orientation.z), 0,
			                  0,                  0, 1
		);
		cameraAngle = cameraAngle * rotMatrix;
	}

	cameraAngle = lookAt(cameraPos, glm::vec3(0, 0, 0));

	translation = glm::vec3(0, 0, 0);
	rotation    = glm::vec3(0, 0, 0);
	orientation = glm::vec3(0, 0, 0);
}

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &translation, glm::vec3 &rotation, glm::vec3 &orientation, glm::vec3 &light, int &mode) {
	float speed = 0.1;
	float angle = 0.5; // radians
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) orientation.y += angle;
		else if (event.key.keysym.sym == SDLK_RIGHT) orientation.y -= angle;
		else if (event.key.keysym.sym == SDLK_UP) orientation.x += angle;
		else if (event.key.keysym.sym == SDLK_DOWN) orientation.x -= angle;
		else if (event.key.keysym.sym == 'j') orientation.z -= angle;
		else if (event.key.keysym.sym == 'l') orientation.z += angle;
		else if (event.key.keysym.sym == 'a') translation.x -= speed;
		else if (event.key.keysym.sym == 'd') translation.x += speed;
		else if (event.key.keysym.sym == 'w') translation.y += speed;
		else if (event.key.keysym.sym == 's') translation.y -= speed;
		else if (event.key.keysym.sym == 'i') translation.z -= speed;
		else if (event.key.keysym.sym == 'k') translation.z += speed;
		else if (event.key.keysym.sym == 'f') rotation.y -= angle;
		else if (event.key.keysym.sym == 'h') rotation.y += angle;
		else if (event.key.keysym.sym == 't') rotation.x -= angle;
		else if (event.key.keysym.sym == 'g') rotation.x += angle;
		else if (event.key.keysym.sym == ';') light.x -= speed;
		else if (event.key.keysym.sym == '#') light.x += speed;
		else if (event.key.keysym.sym == '[') light.y += speed;
		else if (event.key.keysym.sym == '\'') light.y -= speed;
		else if (event.key.keysym.sym == 'z') mode = 0;
		else if (event.key.keysym.sym == 'x') mode = 1;
		else if (event.key.keysym.sym == 'c') mode = 2;
	} else if (event.type == SDL_KEYUP) {
		/*
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) translation.y = 0;
		else if (event.key.keysym.sym == SDLK_DOWN) translation.y = 0;
		else if (event.key.keysym.sym == 'a') translation.x = 0;
		else if (event.key.keysym.sym == 'd') translation.x = 0;
		else if (event.key.keysym.sym == 'w') translation.z = 0;
		else if (event.key.keysym.sym == 's') translation.z = 0;
		*/
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int animate(std::vector<glm::mat3> instructions, std::vector<int> instructionLengths, int tick, glm::vec3 &translation, glm::vec3 &rotation, glm::vec3 &orientation) {
	int shift = 0;
	int i = 0;
	while (i < instructionLengths.size() && shift + instructionLengths[i] <= tick) {
		shift += instructionLengths[i];
		i++;
	}

	if (i == instructionLengths.size()) { i = 0; tick = 0; }
	else tick += 1;

	translation = instructions[i][0];
	rotation    = instructions[i][1];
	orientation = instructions[i][2];

	return tick;
}

int main(int argc, char *argv[]) {

	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	int mode = 0;

	glm::vec3 translation = glm::vec3(0, 0, 0);
	glm::vec3 rotation    = glm::vec3(0, 0, 0);
	glm::vec3 orientation = glm::vec3(0, 0, 0);

	glm::vec3 cameraPos = glm::vec3(0, 0, 3);
	float focalLength = 1;
	float scale = 1000;

	//     right      up forward
	// x       ?       ?       ?
	// y       ?       ?       ?
	// z       ?       ?       ?
	glm::mat3 cameraAngle = glm::mat3(1);

	std::vector<Colour::Colour> colours = readMTL("models/textured-cornell-box.mtl");
	for (int i = 0; i < colours.size(); i++) {
		std::cout << colours.at(i) << std::endl;
	}

	float modelScale = 0.17;
	std::vector<ModelTriangle::ModelTriangle> triangles = readOBJ("models/textured-cornell-box.obj", modelScale, colours);

	glm::vec3 light = glm::vec3(0, 2.5, 0)*modelScale;

	std::vector<glm::vec3> cameraPath;
	cameraPath.push_back(glm::vec3(0, 0, 3));
	cameraPath.push_back(glm::vec3(-3, 0, 3));
	cameraPath.push_back(glm::vec3(-3, 0, -3));
	cameraPath.push_back(glm::vec3(3, 0, -3));
	cameraPath.push_back(glm::vec3(3, 0, 3));

	int pathPos = 0;

	std::vector<glm::mat3> instructions;
	std::vector<int> instructionLengths;

	instructions.push_back(glm::mat3(glm::vec3(0, 0, 0), glm::vec3(0, (2*M_PI)/90, 0), glm::vec3(0, 0, 0)));
	instructionLengths.push_back(90);
	instructions.push_back(glm::mat3(glm::vec3(0, 0, 0), glm::vec3(0, -(2*M_PI)/90, 0), glm::vec3(0, 0, 0)));
	instructionLengths.push_back(90);
	instructions.push_back(glm::mat3(glm::vec3((2*M_PI)/90, 0, 0), glm::vec3(0, 0, (2*M_PI)/90), glm::vec3(0, 0, 0)));
	instructionLengths.push_back(90);
	instructions.push_back(glm::mat3(glm::vec3((2*M_PI)/90, 0, 0), glm::vec3(0, 0, -(2*M_PI)/90), glm::vec3(0, 0, 0)));
	instructionLengths.push_back(90);

	int tick = 0;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, translation, rotation, orientation, light, mode);
		tick = animate(instructions, instructionLengths, tick, translation, rotation, orientation);
		if (tick == 0) mode = (mode+1)%3;
		update(window, cameraPos, cameraAngle, translation, rotation, orientation, cameraPath, pathPos);
		draw(window, mode, cameraPos, cameraAngle, focalLength, scale, triangles, light);
		//window.savePPM("output/"+std::to_string(mode)+"x"+std::to_string(tick)+".ppm");
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
