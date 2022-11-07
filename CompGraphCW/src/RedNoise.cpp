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

#include <unordered_map>

#include<cstring>

#include <stdio.h>

#define WIDTH 500
#define HEIGHT 500

std::vector<double> zBuffer;

class Camera{

	public: 
		glm::vec3 cameraPosition;
		glm::mat3x3 orientation;
		float focalLength;
		float planeScale;

		double angle = 5 * M_PI/180;

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

			std::cout<<this->cameraPosition.z<<" is z value now"<<std::endl;

		}

		void decreaseZ(){

			this->cameraPosition.z -= 1.0;
		    std::cout<<this->cameraPosition.z<<" is z value now"<<std::endl;
		
		}

		void rotateX(){

			glm::mat3x3 rotationMatrix = glm::mat3x3(1.0, 0.0, 0.0, 0.0, cos(angle), -sin(angle), 0.0, sin(angle), cos(angle));
			this->cameraPosition = rotationMatrix * this->cameraPosition;
		}

		void rotateY(){

			glm::mat3x3 rotationMatrix = glm::mat3x3(cos(angle), 0.0, sin(angle), 0.0, 1.0, 0.0, -sin(angle), 0.0, cos(angle));
			this->cameraPosition = rotationMatrix * this->cameraPosition;
		}

		void changeOrientationY(){
			
			glm::mat3x3 rotationMatrix = glm::mat3x3(cos(angle), 0.0, sin(angle), 0.0, 1.0, 0.0, -sin(angle), 0.0, cos(angle));
			this->orientation = rotationMatrix * this->orientation;
		}

		void changeOrientationX(){

			glm::mat3x3 rotationMatrix = glm::mat3x3(1.0, 0.0, 0.0, 0.0, cos(angle), -sin(angle), 0.0, sin(angle), cos(angle));
			this->orientation = rotationMatrix * this->orientation;
		}
};

std::unordered_map<std::string,Colour> readColours(const std::string &filename){

	std::ifstream inputStream(filename);
	std::string line;
	std:: unordered_map<std::string,Colour> colourMap;

	while(inputStream && std::getline(inputStream,line)){

		std::string colourName = line;
		std::vector<std::string> splitColourName = split(colourName,' ');
		std::getline(inputStream,line);
		std::vector<std::string> splitLine = split(line,' ');
		colourMap[splitColourName[1]] =  Colour(stof(splitLine[1])*255, stof(splitLine[2])*255, stof(splitLine[3])*255);
		std::getline(inputStream,line);
	}

	return colourMap;

}

 std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues){

	std::vector<glm::vec3> result;
	glm::vec3 numVal = glm::vec3(numberOfValues-1,numberOfValues-1,numberOfValues-1);
	glm::vec3 stepValues =  (to-from)/numVal;
	result.push_back(from);
	for(int i = 1; i< numberOfValues; i++){

			from = from + stepValues;
			result.push_back(from);
			
		}

		return result;
	
}

CanvasPoint getCanvasIntersectionPoint(Camera &camera, glm::vec3 vertex){

	vertex = -(vertex - camera.cameraPosition);
	vertex = vertex * camera.orientation;

	
	float u = camera.focalLength *  vertex.x/vertex.z;
	float v = camera.focalLength * vertex.y/vertex.z;

	
	u = -u *camera.planeScale+ WIDTH/2.0;
	v = v *camera.planeScale + HEIGHT/2.0;

	CanvasPoint p = CanvasPoint(u,v);
	p.depth = vertex.z;

	return p;
	
	
}

std::vector<ModelTriangle> readOBJFile(const std::string &filename,float scalingFactor,std::unordered_map<std::string,Colour> colourMap,float focalLength,Camera &camera, DrawingWindow &window ){


	std::ifstream inputStream(filename);
	std::string line;

	std::vector<ModelTriangle> objTriangles;
	std::vector<glm::vec3> coordinates;
	Colour currentColour;

	while(std::getline(inputStream,line)){

		std::vector<std::string> splitLine = split(line, ' ');

		if(splitLine[0].find("usemtl") != std::string::npos){
			std::string colourName = splitLine[1];
			currentColour = colourMap[colourName];
		}
		if(splitLine[0] == "v"){
			float x = std::stof(splitLine[1]) * scalingFactor;
			float y = std::stof(splitLine[2]) * scalingFactor;
			float z = std::stof(splitLine[3]) * scalingFactor;
			coordinates.push_back(glm::vec3(x,y,z));
		} 
		
		else if(splitLine[0] == "f"){
			int point1 =std::stoi(splitLine[1])-1;
			int point2 =std::stoi(splitLine[2])-1;
			int point3 =std::stoi(splitLine[3])-1;
			objTriangles.push_back(ModelTriangle(coordinates.at(point1), coordinates.at(point2),coordinates.at(point3),currentColour));


			
		}
	}

	return objTriangles;


}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues){

	float step = (to-from)/(numberOfValues-1);
	std::vector<float> result;
	float generated = from;
	 result.push_back(from);
		for(int i = 1; i< numberOfValues; i++){

			generated = generated + step;
			result.push_back(generated);
			
		}

	return result;

}

 void drawLine(CanvasPoint from, CanvasPoint to, DrawingWindow &window, Colour colour,Camera &camera){	
		
		//float xDiff = std::fmax(to.x,0) - std::fmin(from.x,WIDTH);
		float xDiff = to.x - from.x;
		float yDiff = to.y - from.y;
		float zDiff = to.depth - from.depth;
		float noOfSteps = std::fmax(abs(xDiff), abs(yDiff));
		float xStepSize = xDiff/noOfSteps;
		float yStepSize = yDiff/noOfSteps;
		float zStepSize = zDiff/noOfSteps;
		

		// std::vector<float> zs = interpolateSingleFloats(from.depth,to.depth,noOfSteps);
		// std::vector<float> xs = interpolateSingleFloats(from.x,to.x,noOfSteps);
		// std::vector<float> ys = interpolateSingleFloats(from.y,to.y,noOfSteps);
		float red = colour.red;
		float green = colour.green;
		float blue = colour.blue;
		uint32_t pixelColour = (255 << 24) + (int(red)<<16) + (int(green)<<8 ) + (int(blue));
		
		for( float i = 0.0; i<noOfSteps; i++){

			float x = from.x + (xStepSize * i);
			float y = from.y + (yStepSize * i);
			float z = from.depth + (zStepSize* i);

			
			int ix = int(x);
			int iy = int(y);

			if(ix + iy * WIDTH <=zBuffer.size()){
				
				double prev_depth = zBuffer[ix + iy * WIDTH];
				if ((1.0/z) >= prev_depth && (ix>=0 && ix <WIDTH && iy >=0 && iy <HEIGHT)) {
				window.setPixelColour(ix, iy, pixelColour);
				zBuffer[ix + iy * WIDTH] = (1.0/z);
			}

			}
			
			


		}



		// for(int i = 0; i < noOfSteps;i++){

		// 	float x = xs[i];
		// 	float y = ys[i];
		// 	float z = zs[i];

		// 	int ix = int(x);
		// 	int iy = int(y);

		// 	if(ix + iy * WIDTH <=zBuffer.size()){

		// 		double prev_depth = zBuffer[ix + iy * WIDTH];
		// 		if ((1.0/z) >= prev_depth && (ix>=0 && ix <WIDTH && iy >=0 && iy <HEIGHT)) {
		// 		window.setPixelColour(ix, iy, pixelColour);
		// 		zBuffer[ix + iy * WIDTH] = (1.0/z);
		// 	}
		// 	}
			
		// }



 }

void drawTriangle(CanvasTriangle tr, DrawingWindow &window, Colour colour,Camera &camera){

	drawLine(tr.v0(),tr.v1(),window,colour,camera);
	drawLine(tr.v0(),tr.v2(),window,colour,camera);
	drawLine(tr.v1(),tr.v2(),window,colour,camera);

}

void drawFilledTriangle(CanvasTriangle tr, DrawingWindow &window, Colour colour,Camera &camera){

	CanvasPoint top = tr.v0();
	CanvasPoint bottom = tr.v1();
	CanvasPoint right = tr.v2();
	
	if(top.y > bottom.y){

		std::swap(top, bottom);

	}
	if(top.y > right.y){
		std::swap(top, right);
	}

	if(bottom.y < right.y){
		std::swap(bottom,right);
	}

	float y = right.y;

	float x = top.x - (top.x-bottom.x)*(y-top.y)/(bottom.y-top.y);

	CanvasPoint newPoint = CanvasPoint(x,y);

	drawTriangle(tr,window,colour,camera);
	float  z = top.depth - (top.depth - bottom.depth) * (y-top.y)/(bottom.y-top.y);
	

	for(float Y = std::fmax(0,top.y); Y < std::fmin(y,HEIGHT); Y++){


		float leftX = top.x - (top.x - x)*(Y-top.y)/(y- top.y);
		float rightX = top.x - (top.x - right.x)*(Y - top.y)/(right.y - top.y);

		float leftZ = top.depth - (top.depth - z) * (Y-top.y)/(y- top.y);
		float rightZ = top.depth- (top.depth-right.depth) *(Y - top.y)/(right.y - top.y);

		drawLine(CanvasPoint(leftX,Y,(leftZ)),CanvasPoint(rightX,Y,(rightZ)),window,colour,camera);
	}

	for(float Y = std::fmax(0,newPoint.y); Y <= std::fmin(bottom.y,HEIGHT); Y++){

		float leftX = newPoint.x - (newPoint.x - bottom.x)*(Y-newPoint.y)/(bottom.y-newPoint.y);
		float rightX = right.x - (right.x - bottom.x)*(Y - right.y)/(bottom.y - right.y);

		float leftZ = z - ( z - bottom.depth ) * (Y-newPoint.y)/(bottom.y-newPoint.y);
		float rightZ = right.depth - (right.depth - bottom.depth) *(Y - right.y)/(bottom.y - right.y);

		drawLine(CanvasPoint(leftX,Y, (leftZ)),CanvasPoint(rightX,Y, (rightZ)),window,colour,camera);
	}


}

std::vector<CanvasTriangle> convertModelToCanvas(std::vector<ModelTriangle> vertices, DrawingWindow &window, Camera &camera){

	std::vector<CanvasTriangle> canvasTriangles;

	for(ModelTriangle tr : vertices){


		float z1 = tr.vertices[0].z - camera.cameraPosition.z;
		float z2 = tr.vertices[1].z - camera.cameraPosition.z;
		float z3 = tr.vertices[2].z - camera.cameraPosition.z;

		CanvasPoint p1 = getCanvasIntersectionPoint(camera,tr.vertices[0]);
		CanvasPoint p2 = getCanvasIntersectionPoint(camera,tr.vertices[1]);
		CanvasPoint p3 = getCanvasIntersectionPoint(camera,tr.vertices[2]);		
	
		canvasTriangles.push_back(CanvasTriangle(p1,p2,p3));
	}


	return canvasTriangles;

}



void wireFrameRender(std::vector<ModelTriangle> vertices, DrawingWindow &window, Camera &camera,bool colour) {	

	std::vector<CanvasTriangle> canvasTriangles = convertModelToCanvas(vertices,window,camera);

	window.clearPixels();

	for(int i = 0; i< canvasTriangles.size();i++){
		if(colour==false){

			drawTriangle(canvasTriangles[i],window,Colour(255,255,255),camera);
		} else{

			if(camera.cameraPosition.z != 0)
			drawFilledTriangle(canvasTriangles[i],window,vertices[i].colour,camera);

		}		
	}
}

void drawTexturedTriangle(CanvasTriangle canvasTriangle, CanvasTriangle textureTriangle,DrawingWindow &window, TextureMap textureMap,Camera &camera){


	glm::mat3x3 canvasMatrix = glm::mat3x3(canvasTriangle.v0().x,canvasTriangle.v0().y,1,
										   canvasTriangle.v1().x,canvasTriangle.v1().y,1,
										   canvasTriangle.v2().x,canvasTriangle.v2().y,1);
	
	glm::mat3x3 textureMatrix = glm::mat3x3(textureTriangle.v0().x,textureTriangle.v0().y,1,
											textureTriangle.v1().x,textureTriangle.v1().y,1,
											textureTriangle.v2().x,textureTriangle.v2().y,1);



	glm::mat3x3 mappedTriangle = textureMatrix * (glm::inverse(canvasMatrix));

	CanvasPoint top = canvasTriangle.v0();
	CanvasPoint bottom = canvasTriangle.v1();
	CanvasPoint right = canvasTriangle.v2();
	
	if(top.y > bottom.y){
		std::swap(top, bottom);

	}
	if(top.y > right.y){
		std::swap(top, right);
	}

	if(bottom.y < right.y){
		std::swap(bottom,right);
	}


	
	float x = top.x - (top.x-bottom.x)*(right.y-top.y)/(bottom.y-top.y);

	for(float Y = top.y; Y <right.y; Y++){
	
		float leftX = top.x - (top.x - x)*(Y-top.y)/(right.y- top.y);
		float rightX = top.x - (top.x - right.x)*(Y - top.y)/(right.y - top.y);
		
		if(leftX > rightX){
			std::swap(leftX,rightX);
		}	
		for(float X = leftX; X <rightX; X++){

			glm::vec3 point = mappedTriangle * glm::vec3(X,Y,1);
			window.setPixelColour(X, Y,textureMap.pixels.at(int(point.x) + int(point.y)*textureMap.width));
		}

	}


	for(float Y = right.y; Y <bottom.y; Y++){

		float leftX = x - (x - bottom.x)*(Y-right.y)/(bottom.y-right.y);
		float rightX = right.x - (right.x - bottom.x)*(Y - right.y)/(bottom.y - right.y);

		if(leftX > rightX){
			std::swap(leftX,rightX);
		}

		for(float X = leftX; X <rightX; X++){
			
			glm::vec3 point = mappedTriangle * glm::vec3(X,Y,1);
			
			window.setPixelColour(X, Y,textureMap.pixels.at(int(point.x) + int(point.y)*textureMap.width));
			
		}	
		
		
		}

		drawTriangle(canvasTriangle,window,Colour(255,255,255),camera);	


}

void draw(DrawingWindow &window) {
	window.clearPixels();
	

	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow

	std::vector<glm::vec3> left = interpolateThreeElementValues(topLeft,bottomLeft,HEIGHT);
	std::vector<glm::vec3> right = interpolateThreeElementValues(topRight,bottomRight,HEIGHT);
	for (size_t y = 0; y < window.height; y++) {

		std::vector<glm::vec3> rowInterpol = interpolateThreeElementValues(left.at(y),right.at(y),WIDTH);
		for (size_t x = 0; x < window.width; x++) {
			
			 uint32_t colour = (255 << 24) + (int(rowInterpol.at(x).x) << 16) + ((int(rowInterpol.at(x).y)) << 8) + (int(rowInterpol.at(x).z));
			

			
			 window.setPixelColour(x, y, colour);


		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window,Camera &camera) {
	
	if(event.type == SDL_KEYDOWN){

	 if(event.key.keysym.sym == SDLK_u){

		CanvasPoint p1(rand() % WIDTH, rand() % HEIGHT);
		CanvasPoint p2(rand() % WIDTH, rand() % HEIGHT);
		CanvasPoint p3(rand() % WIDTH, rand() % HEIGHT);
		CanvasTriangle tr(p1,p2,p3);
		Colour colour(rand() % 256, rand() % 256, rand() % 256);
		drawTriangle(tr, window, colour,camera);
	} else if(event.key.keysym.sym == SDLK_f){

		CanvasPoint p1(rand() % WIDTH, rand() % HEIGHT);
		CanvasPoint p2(rand() % WIDTH, rand() % HEIGHT);
		CanvasPoint p3(rand() % WIDTH, rand() % HEIGHT);
		CanvasTriangle tr(p1,p2,p3);
		Colour colour(rand() % 256, rand() % 256, rand() % 256);

		drawFilledTriangle(tr,window,colour,camera);
	} else if(event.key.keysym.sym == SDLK_a){


		std::cout<<camera.cameraPosition.x<<" cameraX initially";
		camera.increaseX();
		std::cout<<" "<<camera.cameraPosition.x<<" and cameraX after"<<std::endl;
	} else if(event.key.keysym.sym == SDLK_d){

		camera.decreaseX();
	} else if(event.key.keysym.sym == SDLK_w){

		camera.decreaseY();
	} else if(event.key.keysym.sym == SDLK_s){

		camera.increaseY();
	} else if(event.key.keysym.sym == SDLK_UP){
		
		camera.decreaseZ();
	} else if(event.key.keysym.sym == SDLK_DOWN){

		camera.increaseZ();
	} else if(event.key.keysym.sym == SDLK_r){

		std::cout<<camera.cameraPosition.x<<" "<<camera.cameraPosition.y<<" "<<camera.cameraPosition.z<<" ->";
		
		std::cout<<camera.cameraPosition.x<<" "<<camera.cameraPosition.y<<" "<<camera.cameraPosition.z<<std::endl;

		camera.rotateX();
		camera.changeOrientationX();
		std::cout<<"That was around the X axis"<<std::endl;



	} else if(event.key.keysym.sym == SDLK_q){

		std::cout<<camera.cameraPosition.x<<" "<<camera.cameraPosition.y<<" "<<camera.cameraPosition.z<<" ->";
		// double angle = 30 * M_PI/180;
		

		// float cosVal = cos(angle);
		// float sinVal = sin(angle);

		// // if(abs(sinVal - 0) < 0.0001){

		// // 	sinVal = 0;
		// // }

		// // if(abs(cosVal - 0) < 0.0001){

		// // 	cosVal = 0 ;
		// // }
		// glm::vec3 col1 = glm::vec3(cosVal,0,-sinVal);
		// glm::vec3 col2 = glm::vec3(0,1,0);
		// glm::vec3 col3 = glm::vec3(sinVal,0,cosVal);

		// glm::mat3x3 rotationMatrix = glm::mat3x3(col1,
		// 										 col2,
		// 										 col3);
		
		// camera.cameraPosition = rotationMatrix * camera.cameraPosition;

		camera.rotateY();
		camera.changeOrientationY();

		std::cout<<camera.cameraPosition.x<<" "<<camera.cameraPosition.y<<" "<<camera.cameraPosition.z<<" "<<std::endl;
		std::cout<<"That was around the y axis"<<std::endl;
	} 
}

}


int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	srand(time(NULL));
	
	zBuffer = std::vector<double>(WIDTH * HEIGHT, 0.0);
	
	glm::vec3 cameraPos = glm::vec3(0.0,0.0,4.0);
	float focalLength = 2.0;
	float planeScale = 250.0;
	Camera camera = Camera(cameraPos, glm::mat3x3(),focalLength,planeScale);


	TextureMap map = TextureMap("src/texture.ppm");
	std::unordered_map<std::string,Colour> colourMap = readColours("src/cornell-box.mtl");
	std::vector<ModelTriangle> modelTriangles = readOBJFile("src/cornell-box.obj",0.35,colourMap,2,camera,window);

	bool colour = true;

	
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window,camera);
		//draw(window);
		wireFrameRender(modelTriangles,window,camera,colour);
		

		zBuffer = std::vector<double>(WIDTH * HEIGHT, 0.0);		

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();

	}

}

