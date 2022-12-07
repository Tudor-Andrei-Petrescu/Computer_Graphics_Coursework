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

#include<unistd.h>   

#define WIDTH 640
#define HEIGHT 480

enum lightType {PROXIMITY,AOI,SPECULAR,AMBIENT};

enum shading {GOURAUD, PHONG,HARD,SOFT};

std::vector<double> zBuffer;




class Camera{

public:
	glm::vec3 cameraPosition;
	glm::mat3x3 orientation;
	glm::vec3 lightSource = glm::vec3(0,0.7,0.5);


	float focalLength;
	float planeScale;

	int lightType = -1;
	int shading = -1;
	
	bool mirror = false;
	bool glass = false;
	bool orbit = false;

		

	double angle = 2 * M_PI / 180;

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

	void increaseLightX(){

		this->lightSource.x += 0.1;
	}

	void decreaseLightX(){

		this->lightSource.x -= 0.1;
	}

	void increaseLightY(){

		this->lightSource.y += 0.1;
	}

	void decreaseLightY(){

		this->lightSource.y -= 0.1;
	}
	
	void resetLight(){
		this->lightSource = glm::vec3(0,0.7,0.5);
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

			tr.vertexIndices[0] = point1+1;
			tr.vertexIndices[1] = point2+1;
			tr.vertexIndices[2] = point3+1;
	
			tr.normal = glm::normalize(glm::cross(tr.vertices[1] - tr.vertices[0], tr.vertices[2] - tr.vertices[0]));
			objTriangles.push_back(tr);
		}

		int size = coordinates.size();

		std::vector<glm::vec3> coordNormals;
		std::vector<int> pointsFrequency;

		for(int i = 0; i <size+2; i++ ){

			coordNormals.push_back(glm::vec3(0,0,0));
			pointsFrequency.push_back(0);
		}

		for (ModelTriangle tr : objTriangles){

			coordNormals[tr.vertexIndices[0]] += tr.normal;
			coordNormals[tr.vertexIndices[1]] += tr.normal;
			coordNormals[tr.vertexIndices[2]] += tr.normal;

			pointsFrequency[tr.vertexIndices[0]]++;
			pointsFrequency[tr.vertexIndices[1]]++;
			pointsFrequency[tr.vertexIndices[2]]++;


			
		}

		for(ModelTriangle tr : objTriangles){

			tr.vertexNormals[0] = glm::normalize(coordNormals[tr.vertexIndices[0]] / (float)pointsFrequency[tr.vertexIndices[0]]);	
			tr.vertexNormals[1] = glm::normalize(coordNormals[tr.vertexIndices[1]] / (float)pointsFrequency[tr.vertexIndices[1]]);
			tr.vertexNormals[2] = glm::normalize(coordNormals[tr.vertexIndices[2]]  / (float)pointsFrequency[tr.vertexIndices[2]]);

			
			
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
	to.x= round(to.x);
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

		if (ix + iy * WIDTH <= zBuffer.size()){

			double prev_depth = zBuffer[ix + iy * WIDTH];
			if ((1.0 / z) >= prev_depth && (ix >= 0 && ix < WIDTH && iy >= 0 && iy < HEIGHT)){
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

	for (float Y = std::fmax(0, top.y); Y < std::fmin(y, HEIGHT); Y++){

		float leftX = top.x - (top.x - x) * (Y - top.y) / (y - top.y);
		float rightX = top.x - (top.x - right.x) * (Y - top.y) / (right.y - top.y);

		float leftZ = top.depth - (top.depth - z) * (Y - top.y) / (y - top.y);
		float rightZ = top.depth - (top.depth - right.depth) * (Y - top.y) / (right.y - top.y);

		drawLine(CanvasPoint(leftX, Y, (leftZ)), CanvasPoint(rightX, Y, (rightZ)), window, colour, camera);
	}

	for (float Y = std::fmax(0, newPoint.y); Y <= std::fmin(bottom.y, HEIGHT); Y++){

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


glm::vec3 getBarycentricCoordinates(glm::vec3 point, ModelTriangle triangle){

	glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
	glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
	glm::vec3 e2 = point - triangle.vertices[0];

	float dot00 = glm::dot(e0, e0);
	float dot01 = glm::dot(e0, e1);
	float dot11 = glm::dot(e1, e1);
	float dot20 = glm::dot(e2, e0);
	float dot21 = glm::dot(e2, e1);
	float denom = dot00 * dot11 - dot01 * dot01;
	float v = (dot11 * dot20 - dot01 * dot21) / denom;
	float w = (dot00 * dot21 - dot01 * dot20) / denom;

	float u = 1.0f - v - w;
	return glm::vec3(u, v, w);

}

float proximity(glm::vec3 point, glm::vec3 light){

	float distance = glm::length(point - light);

	float brightness =(1.0f/(3.0f*M_1_PI*distance*distance));

		if(brightness > 1.0f){

				brightness = 1.0f;
			}

		return brightness;


}

float angleOfIncidence(RayTriangleIntersection closestIntersection, Camera camera){

	float aoi = glm::dot(closestIntersection.intersectedTriangle.normal, glm::normalize(camera.lightSource - closestIntersection.intersectionPoint));


	
	if(aoi < 0.0f || aoi > 1.0f){

		aoi = 1.0f;


	}


	return aoi;

}



float specular(RayTriangleIntersection closestIntersection, Camera camera){

	glm::vec3 reflectedLight = glm::normalize(closestIntersection.intersectionPoint - camera.lightSource) - 2.0f * glm::dot(glm::normalize(closestIntersection.intersectionPoint - camera.lightSource), closestIntersection.intersectedTriangle.normal) * closestIntersection.intersectedTriangle.normal;
	float specular = glm::dot(reflectedLight, glm::normalize(camera.cameraPosition - closestIntersection.intersectionPoint));

				
	specular = pow(specular, 256.0f);

	if(specular < 0.0f){

		specular = 0.0f;
	}

	if(specular > 1.0f){

		specular = 1.0f;
	}

	return specular;

}

std::vector<std::vector<glm::vec3>> getLights(glm::vec3 lightCentre){

	std::vector<std::vector<glm::vec3>> lights;

	for(int i = -2; i < 2; i++){

		std::vector<glm::vec3> light;

		for(int j = -2; j < 2; j++){

			light.push_back(glm::vec3(lightCentre.x + (i*0.05f), lightCentre.y, lightCentre.z+ (j*0.05f)));

		}

		lights.push_back(light);

	}

	return lights;

}



void fastRendering(DrawingWindow &window, Camera &camera,std::vector<ModelTriangle> modelTriangles,CanvasPoint point){

	//LORD PLEASE FORGIVE ME FOR I WHAT I HAVE DONE TO THIS SECTION OF THE CODE 


	glm::vec3 rayDirection = getWorldIntersectionPoint(camera, point);
	RayTriangleIntersection closestIntersection = getClosestIntersectionPoints(camera.cameraPosition, rayDirection, modelTriangles) ;

	if(closestIntersection.triangleIndex != modelTriangles.size()+1){

		

		Colour colour = closestIntersection.intersectedTriangle.colour;

		RayTriangleIntersection shadow = getClosestIntersectionPoints(camera.lightSource,closestIntersection.intersectionPoint - camera.lightSource, modelTriangles);


		float red = int(colour.red);
		float green = int(colour.green);
		float blue = int(colour.blue);
		float brightness = proximity(closestIntersection.intersectionPoint,camera.lightSource);

		if(camera.mirror == true && closestIntersection.intersectedTriangle.colour.name == "newmtl Blue"){
			
			glm::vec3 normal = closestIntersection.intersectedTriangle.normal;
			glm::vec3 rayDirNormal = glm::normalize(closestIntersection.intersectionPoint - camera.cameraPosition);
			glm::vec3 reflected= rayDirNormal - 2.0f * glm::dot(rayDirNormal,normal) * normal;
			RayTriangleIntersection reflection = getClosestIntersectionPoints(closestIntersection.intersectionPoint, reflected, modelTriangles);
			
			if(reflection.intersectedTriangle.colour.name != "newmtl Blue"){
				
			
				red =(int(reflection.intersectedTriangle.colour.red)*0.75) ;
				green = (int(reflection.intersectedTriangle.colour.green)*0.75);
				blue = (int(reflection.intersectedTriangle.colour.blue)*0.75);
			
			}

		}

		

		if(camera.lightType != -1){

		if(camera.lightType == AOI || camera.lightType == SPECULAR){
			float aoi = angleOfIncidence(closestIntersection, camera);
			brightness *= aoi;

			if(camera.lightType == SPECULAR){

			float spec = specular(closestIntersection, camera);
			brightness = brightness + 0.7f*spec;
			
		}
		}
		
		if(camera.lightType == AMBIENT){

		

			if(brightness < 0.3f){

				brightness = 0.3f;
			}
		}

		red *= brightness;
		green *= brightness; 
		blue *= brightness;
		}

		if(camera.shading != -1){

			glm::vec3 barycentric = getBarycentricCoordinates(closestIntersection.intersectionPoint, closestIntersection.intersectedTriangle);

			float u = barycentric.x;
			float v = barycentric.y;
			float w = barycentric.z;


			if(camera.shading == GOURAUD){
				
				glm::vec3 p0 = closestIntersection.intersectedTriangle.vertexNormals[0];
				glm::vec3 p1 = closestIntersection.intersectedTriangle.vertexNormals[1];
				glm::vec3 p2 = closestIntersection.intersectedTriangle.vertexNormals[2];

				glm::vec3 lightNormalised0 = glm::normalize(camera.lightSource-closestIntersection.intersectedTriangle.vertices[0]);
				glm:: vec3 lightNormalised1 = glm::normalize(camera.lightSource-closestIntersection.intersectedTriangle.vertices[1]);
				glm::vec3 lightNormalised2 = glm::normalize(camera.lightSource-closestIntersection.intersectedTriangle.vertices[2]);


				glm::vec3 cameraNormalised0 = glm::normalize(camera.cameraPosition - closestIntersection.intersectedTriangle.vertices[0]);
				glm::vec3 cameraNormalised1 = glm::normalize(camera.cameraPosition - closestIntersection.intersectedTriangle.vertices[1]);
				glm::vec3 cameraNormalised2 = glm::normalize(camera.cameraPosition - closestIntersection.intersectedTriangle.vertices[2]);


				float aoi0 = glm::dot(closestIntersection.intersectedTriangle.vertices[0],lightNormalised0);
				float aoi1 = glm::dot(closestIntersection.intersectedTriangle.vertices[1],lightNormalised1);
				float aoi2 = glm::dot(closestIntersection.intersectedTriangle.vertices[2],lightNormalised2);


				if(aoi0 <0.0f || aoi0 > 1.0f){
					aoi0 = 1.0f;

				}
				if(aoi1 <0.0f || aoi1 > 1.0f){
					aoi1 = 1.0f;

				}
				if(aoi2 <0.0f || aoi2 > 1.0f){
					aoi2 = 1.0f;

				}


				glm::vec3 reflectedLight0 = glm::normalize(lightNormalised0 - 2.0f * glm::dot(lightNormalised0, p0) * p0);
				glm::vec3 reflectedLight1 = glm::normalize(lightNormalised1 - 2.0f * glm::dot(lightNormalised1, p1) * p1);
				glm::vec3 reflectedLight2 = glm::normalize(lightNormalised2 - 2.0f * glm::dot(lightNormalised2, p2) * p2);

				float spec0 = glm::dot(reflectedLight0, cameraNormalised0);
			    float spec1 = glm::dot(reflectedLight1, cameraNormalised1);
				float spec2 = glm::dot(reflectedLight2, cameraNormalised2);

				spec0 = pow(spec0, 256.0f);
				spec1 = pow(spec1, 256.0f);
				spec2 = pow(spec2, 256.0f);

				if(spec0 <0.0f){
					spec0 = 0.0f;
				}
				if(spec1 <0.0f){
					spec1 = 0.0f;
				}
				if(spec2 <0.0f){
					spec2 = 0.0f;
				}
				if(spec0>1.0f){
					spec0 = 1.0f;
				}
				if(spec1>1.0f){
					spec1 = 1.0f;
				}
				if(spec2>1.0f){
					spec2 = 1.0f;
				}


				float prox0 = proximity(closestIntersection.intersectedTriangle.vertices[0],camera.lightSource);
				float prox1 = proximity(closestIntersection.intersectedTriangle.vertices[1],camera.lightSource);
				float prox2 = proximity(closestIntersection.intersectedTriangle.vertices[2],camera.lightSource);

				float gouraud =  u * (aoi0 + prox0+spec0) + v * (aoi1 +prox1 +spec1) + w* (aoi2 +prox2 + spec2);
				
				if(gouraud > 1.0f){

					gouraud = 1.0f;
				}
				
				red*=gouraud;
				green*=gouraud;
				blue*=gouraud;

			}
			else if(camera.shading ==PHONG){

			ModelTriangle tr = closestIntersection.intersectedTriangle;
			glm::vec3 pointNormal = tr.vertexNormals[0] * u + tr.vertexNormals[1] * v + tr.vertexNormals[2] * w;


			float aoi = glm::dot(closestIntersection.intersectionPoint, glm::normalize(camera.lightSource - closestIntersection.intersectionPoint));

			if(aoi < 0.0f || aoi > 1.0f){

				aoi = 1.0f;

			}

			glm::vec3 reflectedLight = glm::normalize(closestIntersection.intersectionPoint - camera.lightSource) - 2.0f * glm::dot(glm::normalize(closestIntersection.intersectionPoint - camera.lightSource), pointNormal) * pointNormal;
			float specular = glm::dot(reflectedLight, glm::normalize(camera.cameraPosition - closestIntersection.intersectionPoint));

				
			specular = pow(specular, 256.0f);

			if(specular < 0.0f){

				specular = 0.0f;
			}

			if(specular > 1.0f){

				specular = 1.0f;
			}	

			float prox = proximity(closestIntersection.intersectionPoint, camera.lightSource);

			float phong = aoi + prox *specular;

			if(phong > 1.0f){

				phong = 1.0f;
			}


			red *= phong;
			green *= phong;
			blue *=phong;

			

		
		}
		else if (camera.shading == HARD){

			if(shadow.triangleIndex != modelTriangles.size()+1 && shadow.triangleIndex != closestIntersection.triangleIndex){


				red*=0.6f;
				green*=0.6f;
				blue*=0.6f;
			}
		}
		else if(camera.shading == SOFT){

			std::vector<std::vector<glm::vec3>> lightPoints = getLights(camera.lightSource);


			for(int i = 0; i < lightPoints.size(); i++){

				for(int j = 0; j < lightPoints[i].size();j++){
					float intensity = 1.0f;
					RayTriangleIntersection shadowRay = getClosestIntersectionPoints(lightPoints.at(i).at(j),closestIntersection.intersectionPoint - lightPoints.at(i).at(j), modelTriangles);

					if(shadowRay.triangleIndex != modelTriangles.size()+1 && shadowRay.triangleIndex != closestIntersection.triangleIndex){

						intensity-=0.04f;
					}

					red*=intensity;
					green*=intensity;
					blue*=intensity;

				}
			}
		}

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


void drawRayTraced(DrawingWindow &window, Camera &camera, std::vector<ModelTriangle> modelTriangles){

	window.clearPixels();

	std::vector<std::thread> threads = std::vector<std::thread>();

	for(int i = 0; i<4;i++){

			std::thread t(computePixelValues, i*120, (i+1)*120, std::ref(window), std::ref(camera), modelTriangles);
			threads.push_back(std::move(t));
			
		}

	for(auto &t : threads){

			t.join();
		}

	

}


void draw(std::vector<ModelTriangle> vertices, DrawingWindow &window, Camera &camera, bool colour,bool lookAt, bool raytrace,bool sphere,std::vector<ModelTriangle> sphereTr)
{

	std::vector<CanvasTriangle> canvasTriangles = convertModelToCanvas(vertices, window, camera);

	window.clearPixels();

	if(lookAt == true || camera.orbit == true){

		camera.rotateY(true);

		if(lookAt == true)
		camera.lookAt(glm::vec3(0,0,0));
	}
	


	if(raytrace == true){

		if(sphere){
			drawRayTraced(window, camera, sphereTr);
			return;
		}

		drawRayTraced(window, camera,vertices);
		return;
	}

	for (int i = 0; i < canvasTriangles.size(); i++){

		if (colour == false){

			drawTriangle(canvasTriangles[i], window, Colour(255, 255, 255), camera);
		}
		else{

			if (camera.cameraPosition.z !=0)

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

void handleEvent(SDL_Event event, DrawingWindow &window, Camera &camera,bool &colour, bool &lookAt, bool &raytrace,bool &sphere)
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

			lookAt = true;
		}
		else if (event.key.keysym.sym == SDLK_u){
			
			lookAt = false;
			camera.orbit = false;

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
		else if(event.key.keysym.sym == SDLK_5){

			if(camera.lightType == SPECULAR){

				camera.lightType = -1;
				std::cout<<"SPECULAR light off"<<std::endl;
			}

			else{
				camera.lightType = SPECULAR;
				std::cout<<"SPECULAR light on"<<std::endl;
			}
		}
		else if(event.key.keysym.sym == SDLK_6){

			if(camera.lightType == AMBIENT){

				camera.lightType = -1;
				std::cout<<"AMBIENT light off"<<std::endl;
			}

			else{
				camera.lightType = AMBIENT;
				std::cout<<"AMBIENT light on"<<std::endl;
			}
		}
		else if(event.key.keysym.sym == SDLK_p){

			std::cout<<"Sphere on"<<std::endl;

			sphere = !sphere;
		}
		else if(event.key.keysym.sym == SDLK_g){
			

			if(camera.shading == -1){
			std::cout<<"Gourad on"<<std::endl;

			camera.shading = GOURAUD;
			camera.lightSource = glm::vec3(-1,3.1,1);
			}
			else{
				std::cout<<"Gourad off"<<std::endl;
				camera.shading = -1;
				camera.lightSource = glm::vec3(0,0.7,0.5);
			}

		}
		else if(event.key.keysym.sym ==SDLK_7){

				if(camera.shading == -1){
			std::cout<<"Phong on"<<std::endl;

			camera.shading = PHONG;
			camera.lightSource = glm::vec3(-1,3.1,1);
			}
			else{
				std::cout<<"Phong off"<<std::endl;
				camera.shading = -1;
				camera.lightSource = glm::vec3(0,0.7,0.5);
			}

		}
		else if(event.key.keysym.sym == SDLK_LSHIFT){
			camera.increaseLightX();
			std::cout<<"Light is now at:"<<camera.lightSource.x<<","<<camera.lightSource.y<<","<<camera.lightSource.z<<std::endl;
		}
		else if(event.key.keysym.sym == SDLK_RSHIFT){
			camera.decreaseLightX();
			std::cout<<"Light is now at:"<<camera.lightSource.x<<","<<camera.lightSource.y<<","<<camera.lightSource.z<<std::endl;
		}
		else if(event.key.keysym.sym == SDLK_LCTRL){
			camera.increaseLightY();
			std::cout<<"Light is now at:"<<camera.lightSource.x<<","<<camera.lightSource.y<<","<<camera.lightSource.z<<std::endl;
		}
		else if(event.key.keysym.sym == SDLK_RCTRL){
			camera.decreaseLightY();
			std::cout<<"Light is now at:"<<camera.lightSource.x<<","<<camera.lightSource.y<<","<<camera.lightSource.z<<std::endl;
		}
		else if(event.key.keysym.sym == SDLK_z){
			
			camera.mirror = !camera.mirror;

			
		}
		else if(event.key.keysym.sym == SDLK_x){

			camera.glass = !camera.glass;
		}
		else if(event.key.keysym.sym == SDLK_8){

			if(camera.shading == -1){
			std::cout<<"HARD SHADING on"<<std::endl;

			camera.shading = HARD;
			}
			else{
				std::cout<<"HARD SHADING off"<<std::endl;
				camera.shading = -1;
				
			}
		}
		else if(event.key.keysym.sym == SDLK_c){

				if(camera.shading == -1){
			std::cout<<"SOFT SHADING on"<<std::endl;

			camera.shading = SOFT;
			}
			else{
				std::cout<<"SOFT SHADING off"<<std::endl;
				camera.shading = -1;
				
			}	

		}
		else if(event.key.keysym.sym == SDLK_9){

			camera.orbit = true;
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
	
	std::vector<ModelTriangle> sphereTr = readOBJFile("src/sphere.obj", 0.35, colourMap, 2, camera, window);

	bool colour = true;
	bool lookAt = false;
	bool raytrace = false;

	bool sphere = false;

		
	while (true){

		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event))
			handleEvent(event, window, camera,colour, lookAt,raytrace,sphere);
	
		draw(modelTriangles, window, camera, colour, lookAt,raytrace,sphere,sphereTr);
		
		//drawRayTraced(window, camera, modelTriangles);
		zBuffer = std::vector<double>(WIDTH * HEIGHT, 0.0);

		

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();

	
	}
}

