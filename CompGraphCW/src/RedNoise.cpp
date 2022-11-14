#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#include <CanvasPoint.h>
#include <TextureMap.h>
#include <Colour.h>
#include <ModelTriangle.h>

#include <time.h>

#include<limits.h>

#include <unordered_map>

#include <cstring>

#include <stdio.h>

#include <RayTriangleIntersection.h>

#include<thread>

#define WIDTH 500
#define HEIGHT 500

enum lightType {PROXIMITY,AOI,SPECULAR,AMBIENT,GOURAUD};

std::vector<double> zBuffer;




class Camera{

public:
	glm::vec3 cameraPosition;
	glm::mat3x3 orientation;
	glm::vec3 lightSource = glm::vec3(0,0.7,0);


	float focalLength;
	float planeScale;

	int lightType = -1;
	

	double angle = 5 * M_PI / 180;

	Camera(glm::vec3 cameraPosition, glm::mat3x3 orientation, float focalLength, float planeScale){

		this->cameraPosition = cameraPosition;
		this->orientation = orientation;
		this->focalLength = focalLength;
		this->planeScale = planeScale;
	}

	void increaseX(){

		this->cameraPosition.x += 1.0;
	}

	void decreaseX(){

		this->cameraPosition.x -= 1.0;
	}

	void increaseY(){

		this->cameraPosition.y += 1.0;
	}

	void decreaseY(){

		this->cameraPosition.y -= 1.0;
	}

	void increaseZ(){

		this->cameraPosition.z += 1.0;

	}

	void decreaseZ(){

		this->cameraPosition.z -= 1.0;
	}

	void rotateX(bool clockWise){

		double dir = angle;
		if (!clockWise)
			dir = -angle;

		glm::mat3x3 rotationMatrix = glm::mat3x3(1.0, 0.0, 0.0, 0.0, cos(dir), -sin(dir), 0.0, sin(dir), cos(dir));
		this->cameraPosition = rotationMatrix * this->cameraPosition;
	}

	void rotateY(bool clockWise){

		double dir = angle;
		if (!clockWise)
			dir = -angle;
		glm::mat3x3 rotationMatrix = glm::mat3x3(cos(dir), 0.0, sin(dir), 0.0, 1.0, 0.0, -sin(dir), 0.0, cos(dir));
		this->cameraPosition = rotationMatrix * this->cameraPosition;
	}

	void changeOrientationY(bool clockWise){
		double dir = angle;
		if(!clockWise) dir = -angle;
		glm::mat3x3 rotationMatrix = glm::mat3x3(cos(dir), 0.0, sin(dir), 0.0, 1.0, 0.0, -sin(dir), 0.0, cos(dir));
		this->orientation = rotationMatrix * this->orientation;
	}

	void changeOrientationX(bool clockWise){
		double dir = angle;
		if(!clockWise) dir = -angle;
		glm::mat3x3 rotationMatrix = glm::mat3x3(1.0, 0.0, 0.0, 0.0, cos(dir), -sin(dir), 0.0, sin(dir), cos(dir));
		this->orientation = rotationMatrix * this->orientation;
	}

	void lookAt(glm::vec3 point){

		glm::vec3 z = glm::normalize(point-cameraPosition);
		glm::vec3 x = glm::normalize(glm::cross(z,glm::vec3(0,1,0)));
		glm::vec3 y = glm::cross(x,z);
		orientation = glm::mat3x3(x,y,-z);
	}
};

std::unordered_map<std::string, Colour> readColours(const std::string &filename){

	std::ifstream inputStream(filename);
	std::string line;
	std::unordered_map<std::string, Colour> colourMap;

	while (inputStream && std::getline(inputStream, line)){

		std::string colourName = line;
		std::vector<std::string> splitColourName = split(colourName, ' ');
		std::getline(inputStream, line);
		std::vector<std::string> splitLine = split(line, ' ');

		Colour colour = Colour(stof(splitLine[1]) * 255, stof(splitLine[2]) * 255, stof(splitLine[3]) * 255);
		colour.name = colourName;
		colourMap[splitColourName[1]] = colour;
		std::getline(inputStream, line);
	}

	return colourMap;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues)
{

	std::vector<glm::vec3> result;
	glm::vec3 numVal = glm::vec3(numberOfValues - 1, numberOfValues - 1, numberOfValues - 1);
	glm::vec3 stepValues = (to - from) / numVal;
	result.push_back(from);
	for (int i = 1; i < numberOfValues; i++){

		from = from + stepValues;
		result.push_back(from);
	}

	return result;
}

CanvasPoint getCanvasIntersectionPoint(Camera &camera, glm::vec3 vertex)
{

	vertex = -(vertex - camera.cameraPosition);
	vertex = vertex * camera.orientation;

	float u = camera.focalLength * vertex.x / vertex.z;
	float v = camera.focalLength * vertex.y / vertex.z;


	u = -u * camera.planeScale + WIDTH / 2.0;
	v = v * camera.planeScale + HEIGHT / 2.0;

	CanvasPoint p = CanvasPoint(u, v);
	p.depth = vertex.z;


	return p;
}

glm::vec3 getWorldIntersectionPoint(Camera &camera, CanvasPoint point)
{
	
	float u = point.x - WIDTH/2.0;
	float v = HEIGHT/2.0- point.y;
	
	float z = -camera.focalLength*camera.planeScale;

	glm::vec3 worldPoint = camera.orientation * glm::vec3(u, v, z);


	return glm::normalize(worldPoint);

	
}

std::vector<ModelTriangle> readOBJFile(const std::string &filename, float scalingFactor, std::unordered_map<std::string, Colour> colourMap, float focalLength, Camera &camera, DrawingWindow &window){

	std::ifstream inputStream(filename);
	std::string line;

	std::vector<ModelTriangle> objTriangles;
	std::vector<glm::vec3> coordinates;
	Colour currentColour;

	while (std::getline(inputStream, line)){

		std::vector<std::string> splitLine = split(line, ' ');

		if (splitLine[0].find("usemtl") != std::string::npos){

			std::string colourName = splitLine[1];
			currentColour = colourMap[colourName];
		}
		if (splitLine[0] == "v"){

			float x = std::stof(splitLine[1]) * scalingFactor;
			float y = std::stof(splitLine[2]) * scalingFactor;
			float z = std::stof(splitLine[3]) * scalingFactor;
			coordinates.push_back(glm::vec3(x, y, z));
		}
		else if (splitLine[0] == "f"){
			
			int point1 = std::stoi(splitLine[1]) - 1;
			int point2 = std::stoi(splitLine[2]) - 1;
			int point3 = std::stoi(splitLine[3]) - 1;

			ModelTriangle tr = ModelTriangle(coordinates.at(point1), coordinates.at(point2), coordinates.at(point3), currentColour);

			tr.normal = glm::normalize(glm::cross(tr.vertices[1] - tr.vertices[0], tr.vertices[2] - tr.vertices[0]));
			objTriangles.push_back(tr);
		}
	}

	return objTriangles;
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues){

	float step = (to - from) / (numberOfValues - 1);
	std::vector<float> result;
	float generated = from;
	result.push_back(from);
	for (int i = 1; i < numberOfValues; i++){

		generated = generated + step;
		result.push_back(generated);
	}

	return result;
}

void drawLine(CanvasPoint from, CanvasPoint to, DrawingWindow &window, Colour colour, Camera &camera){

	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float zDiff = to.depth - from.depth;

	float noOfSteps = std::fmax(abs(xDiff), abs(yDiff));
	
	float xStepSize = xDiff / noOfSteps;
	float yStepSize = yDiff / noOfSteps;
	float zStepSize = zDiff / noOfSteps;

	float red = colour.red;
	float green = colour.green;
	float blue = colour.blue;
	uint32_t pixelColour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + (int(blue));

	for (float i = 0.0; i < noOfSteps; i++){

		float x = from.x + (xStepSize * i);
		float y = from.y + (yStepSize * i);
		float z = from.depth + (zStepSize * i);

		int ix = int(x);
		int iy = int(y);

		if (ix + iy * WIDTH <= zBuffer.size())
		{

			double prev_depth = zBuffer[ix + iy * WIDTH];
			if ((1.0 / z) >= prev_depth && (ix >= 0 && ix < WIDTH && iy >= 0 && iy < HEIGHT))
			{
				window.setPixelColour(ix, iy, pixelColour);
				zBuffer[ix + iy * WIDTH] = (1.0 / z);
			}
		}
	}


}

void drawTriangle(CanvasTriangle tr, DrawingWindow &window, Colour colour, Camera &camera){

	drawLine(tr.v0(), tr.v1(), window, colour, camera);
	drawLine(tr.v0(), tr.v2(), window, colour, camera);
	drawLine(tr.v1(), tr.v2(), window, colour, camera);
}

void drawFilledTriangle(CanvasTriangle tr, DrawingWindow &window, Colour colour, Camera &camera){

	CanvasPoint top = tr.v0();
	CanvasPoint bottom = tr.v1();
	CanvasPoint right = tr.v2();

	if (top.y > bottom.y){

		std::swap(top, bottom);
	}
	if (top.y > right.y){

		std::swap(top, right);
	}
	if (bottom.y < right.y){

		std::swap(bottom, right);
	}

	float y = right.y;

	float x = top.x - (top.x - bottom.x) * (y - top.y) / (bottom.y - top.y);

	CanvasPoint newPoint = CanvasPoint(x, y);

	drawTriangle(tr, window, colour, camera);

	float z = top.depth - (top.depth - bottom.depth) * (y - top.y) / (bottom.y - top.y);

	for (float Y = std::fmax(0, top.y); Y < std::fmin(y, HEIGHT); Y++)
	{

		float leftX = top.x - (top.x - x) * (Y - top.y) / (y - top.y);
		float rightX = top.x - (top.x - right.x) * (Y - top.y) / (right.y - top.y);

		float leftZ = top.depth - (top.depth - z) * (Y - top.y) / (y - top.y);
		float rightZ = top.depth - (top.depth - right.depth) * (Y - top.y) / (right.y - top.y);

		drawLine(CanvasPoint(leftX, Y, (leftZ)), CanvasPoint(rightX, Y, (rightZ)), window, colour, camera);
	}

	for (float Y = std::fmax(0, newPoint.y); Y <= std::fmin(bottom.y, HEIGHT); Y++)
	{

		float leftX = newPoint.x - (newPoint.x - bottom.x) * (Y - newPoint.y) / (bottom.y - newPoint.y);
		float rightX = right.x - (right.x - bottom.x) * (Y - right.y) / (bottom.y - right.y);

		float leftZ = z - (z - bottom.depth) * (Y - newPoint.y) / (bottom.y - newPoint.y);
		float rightZ = right.depth - (right.depth - bottom.depth) * (Y - right.y) / (bottom.y - right.y);

		drawLine(CanvasPoint(leftX, Y, (leftZ)), CanvasPoint(rightX, Y, (rightZ)), window, colour, camera);
	}
}

std::vector<CanvasTriangle> convertModelToCanvas(std::vector<ModelTriangle> vertices, DrawingWindow &window, Camera &camera)
{

	std::vector<CanvasTriangle> canvasTriangles;

	for (ModelTriangle tr : vertices){

		CanvasPoint p1 = getCanvasIntersectionPoint(camera, tr.vertices[0]);
		CanvasPoint p2 = getCanvasIntersectionPoint(camera, tr.vertices[1]);
		CanvasPoint p3 = getCanvasIntersectionPoint(camera, tr.vertices[2]);

		canvasTriangles.push_back(CanvasTriangle(p1, p2, p3));
	}

	return canvasTriangles;
}


RayTriangleIntersection getClosestIntersectionPoints(glm::vec3 cameraPosition, glm::vec3 rayDirection,std::vector<ModelTriangle> modelTriangles){

	RayTriangleIntersection closestIntersection;

	closestIntersection.triangleIndex = modelTriangles.size()+1;
	float minDist = 999999999.0f;

	for(int i = 0; i < modelTriangles.size(); i++){

		glm::vec3 e0 = modelTriangles[i].vertices[1] - modelTriangles[i].vertices[0];
		glm::vec3 e1 = modelTriangles[i].vertices[2] - modelTriangles[i].vertices[0];
		glm::vec3 SPVector = cameraPosition - modelTriangles[i].vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

		float t = possibleSolution.x;
		float u = possibleSolution.y;
		float v = possibleSolution.z;
		
	
		if(t < minDist  && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) &&(u + v) <= 1.0  && t>=0.001f){

			minDist = t;

			closestIntersection.intersectionPoint = cameraPosition + (t * rayDirection);
			closestIntersection.intersectedTriangle = modelTriangles[i];
			closestIntersection.distanceFromCamera = t;
			closestIntersection.triangleIndex = i;




		}


	}

	return closestIntersection;



}

void fastRendering(DrawingWindow &window, Camera &camera,std::vector<ModelTriangle> modelTriangles,CanvasPoint point){


	glm::vec3 rayDirection = getWorldIntersectionPoint(camera, point);
	RayTriangleIntersection closestIntersection = getClosestIntersectionPoints(camera.cameraPosition, rayDirection, modelTriangles) ;

	if(closestIntersection.triangleIndex != modelTriangles.size()+1){

		int index = closestIntersection.triangleIndex;
		Colour colour = closestIntersection.intersectedTriangle.colour;

		RayTriangleIntersection shadow = getClosestIntersectionPoints(camera.lightSource,closestIntersection.intersectionPoint - camera.lightSource, modelTriangles);

		float distance = glm::length(closestIntersection.intersectionPoint - camera.lightSource);

		float red = int(colour.red);
		float green = int(colour.green);
		float blue = int(colour.blue);
		float brightness = (2.0f/(7.0f*M_1_PI*distance*distance));

		if(camera.lightType == PROXIMITY || camera.lightType == AOI){

			
			if(brightness > 1.0f){

				brightness = 1.0f;
			}

			if(camera.lightType == AOI){

				float aoi = glm::dot(closestIntersection.intersectedTriangle.normal, glm::normalize(camera.lightSource - closestIntersection.intersectionPoint));
			
				if(aoi < 0.0f || aoi > 1.0f){

					aoi = 1.0f;

				}

			brightness *= aoi;

			}
			
			red *= brightness;
			green *= brightness; 
			blue *= brightness;
		}

	

		uint32_t pixelColour = (255 << 24) + ( int(red) << 16)+ (int(green) << 8) + (int(blue));
		window.setPixelColour(point.x,point.y,pixelColour);

	}
}


void computePixelValues(int from, int to, DrawingWindow &window, Camera &camera, std::vector<ModelTriangle> modelTriangles){

	for(int i = from; i < to; i++){
		for(int j = 0; j < window.width; j++){
			CanvasPoint point(j,i);
			fastRendering(window, camera, modelTriangles, point);
		}
	}
}


void drawRasterisedScene(DrawingWindow &window, Camera &camera, std::vector<ModelTriangle> modelTriangles){

	window.clearPixels();

	std::vector<std::thread> threads = std::vector<std::thread>();

	for(int i = 0; i<4;i++){

			std::thread t(computePixelValues, i*125, (i+1)*125, std::ref(window), std::ref(camera), modelTriangles);
			threads.push_back(std::move(t));
			
		}

	for(auto &t : threads){

			t.join();
		}

	

}


void draw(std::vector<ModelTriangle> vertices, DrawingWindow &window, Camera &camera, bool colour,bool orbit, bool raytrace)
{

	std::vector<CanvasTriangle> canvasTriangles = convertModelToCanvas(vertices, window, camera);

	window.clearPixels();

	if(orbit == true){

		camera.rotateY(true);

		camera.lookAt(glm::vec3(0,0,0));
	}

	if(raytrace == true){

		drawRasterisedScene(window, camera,vertices);
		return;
	}

	for (int i = 0; i < canvasTriangles.size(); i++){

		if (colour == false){

			drawTriangle(canvasTriangles[i], window, Colour(255, 255, 255), camera);
		}
		else{

			if (camera.cameraPosition.z != 0)

				drawFilledTriangle(canvasTriangles[i], window, vertices[i].colour, camera);
		}
	}
}

void drawTexturedTriangle(CanvasTriangle canvasTriangle, CanvasTriangle textureTriangle, DrawingWindow &window, TextureMap textureMap, Camera &camera){

	glm::mat3x3 canvasMatrix = glm::mat3x3(canvasTriangle.v0().x, canvasTriangle.v0().y, 1,
										   canvasTriangle.v1().x, canvasTriangle.v1().y, 1,
										   canvasTriangle.v2().x, canvasTriangle.v2().y, 1);

	glm::mat3x3 textureMatrix = glm::mat3x3(textureTriangle.v0().x, textureTriangle.v0().y, 1,
											textureTriangle.v1().x, textureTriangle.v1().y, 1,
											textureTriangle.v2().x, textureTriangle.v2().y, 1);

	glm::mat3x3 mappedTriangle = textureMatrix * (glm::inverse(canvasMatrix));

	CanvasPoint top = canvasTriangle.v0();
	CanvasPoint bottom = canvasTriangle.v1();
	CanvasPoint right = canvasTriangle.v2();

	if (top.y > bottom.y){

		std::swap(top, bottom);
	}
	if (top.y > right.y){

		std::swap(top, right);
	}
	if (bottom.y < right.y){

		std::swap(bottom, right);
	}

	float x = top.x - (top.x - bottom.x) * (right.y - top.y) / (bottom.y - top.y);

	for (float Y = top.y; Y < right.y; Y++){

		float leftX = top.x - (top.x - x) * (Y - top.y) / (right.y - top.y);
		float rightX = top.x - (top.x - right.x) * (Y - top.y) / (right.y - top.y);

		if (leftX > rightX){

			std::swap(leftX, rightX);
		}

		for (float X = leftX; X < rightX; X++){

			glm::vec3 point = mappedTriangle * glm::vec3(X, Y, 1);
			window.setPixelColour(X, Y, textureMap.pixels.at(int(point.x) + int(point.y) * textureMap.width));
		}
	}

	for (float Y = right.y; Y < bottom.y; Y++){

		float leftX = x - (x - bottom.x) * (Y - right.y) / (bottom.y - right.y);
		float rightX = right.x - (right.x - bottom.x) * (Y - right.y) / (bottom.y - right.y);

		if (leftX > rightX){

			std::swap(leftX, rightX);
		}

		for (float X = leftX; X < rightX; X++){

			glm::vec3 point = mappedTriangle * glm::vec3(X, Y, 1);

			window.setPixelColour(X, Y, textureMap.pixels.at(int(point.x) + int(point.y) * textureMap.width));
		}
	}

	drawTriangle(canvasTriangle, window, Colour(255, 255, 255), camera);
}

void handleEvent(SDL_Event event, DrawingWindow &window, Camera &camera,bool &colour, bool &orbit, bool &raytrace)
{

	if (event.type == SDL_KEYDOWN){

		if (event.key.keysym.sym == SDLK_0){

			CanvasPoint p1(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint p2(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint p3(rand() % WIDTH, rand() % HEIGHT);
			CanvasTriangle tr(p1, p2, p3);
			Colour colour(rand() % 256, rand() % 256, rand() % 256);
			drawTriangle(tr, window, colour, camera);
		}
		else if (event.key.keysym.sym == SDLK_f){

			CanvasPoint p1(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint p2(rand() % WIDTH, rand() % HEIGHT);
			CanvasPoint p3(rand() % WIDTH, rand() % HEIGHT);

			CanvasTriangle tr(p1, p2, p3);

			Colour colour(rand() % 256, rand() % 256, rand() % 256);

			drawFilledTriangle(tr, window, colour, camera);
		}
		else if (event.key.keysym.sym == SDLK_a){

			camera.increaseX();
		}
		else if (event.key.keysym.sym == SDLK_d){

			camera.decreaseX();
		}
		else if (event.key.keysym.sym == SDLK_w){

			camera.decreaseY();
		}
		else if (event.key.keysym.sym == SDLK_s){

			camera.increaseY();
		}
		else if (event.key.keysym.sym == SDLK_UP){

			camera.decreaseZ();
		}
		else if (event.key.keysym.sym == SDLK_DOWN){

			camera.increaseZ();
		}
		else if (event.key.keysym.sym == SDLK_k){

			camera.rotateX(true);
			
			camera.lookAt(glm::vec3(0,0,0));
		
		}
		else if (event.key.keysym.sym == SDLK_o){

			orbit = true;
		}
		else if (event.key.keysym.sym == SDLK_u){
			
			orbit = false;

			camera.orientation = glm::mat3x3();

			camera.cameraPosition = glm::vec3(0,0,4);
		}
		else if (event.key.keysym.sym == SDLK_i){

			camera.rotateX(false);

			camera.lookAt(glm::vec3(0,0,0));
		}
		else if (event.key.keysym.sym == SDLK_1){

			colour = !colour;
		}
		else if(event.key.keysym.sym == SDLK_2){

			raytrace = !raytrace;

			std::cout<<"Raytrace on"<<std::endl;
		}
		else if(event.key.keysym.sym == SDLK_3){


			if(camera.lightType == PROXIMITY){

				camera.lightType = -1;
				std::cout<<"Proximity light off"<<std::endl;
			}

			else{
				camera.lightType = PROXIMITY;
				std::cout<<"Proximity light on"<<std::endl;
			}
		}
		else if(event.key.keysym.sym == SDLK_4){

			if(camera.lightType == AOI){

				camera.lightType = -1;
				std::cout<<"AOI light off"<<std::endl;
			}

			else{
				camera.lightType = AOI;
				std::cout<<"AOI light on"<<std::endl;
			}
		}
	}
}


int main(int argc, char *argv[]){

	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	srand(time(NULL));

	zBuffer = std::vector<double>(WIDTH * HEIGHT, 0.0);

	glm::vec3 cameraPos = glm::vec3(0.0, 0.0, 4.0);
	float focalLength = 2.0;
	float planeScale = 250.0;
	Camera camera = Camera(cameraPos, glm::mat3x3(), focalLength, planeScale);

	TextureMap map = TextureMap("src/texture.ppm");
	std::unordered_map<std::string, Colour> colourMap = readColours("src/cornell-box.mtl");
	std::vector<ModelTriangle> modelTriangles = readOBJFile("src/cornell-box.obj", 0.35, colourMap, 2, camera, window);

	bool colour = true;
	bool orbit = false;
	bool raytrace = false;
	
	while (true){

		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event))
			handleEvent(event, window, camera,colour, orbit,raytrace);
	
		draw(modelTriangles, window, camera, colour, orbit,raytrace);
		
		//drawRasterisedScene(window, camera, modelTriangles);
		zBuffer = std::vector<double>(WIDTH * HEIGHT, 0.0);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}

