#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#include <CanvasPoint.h>
#include <TextureMap.h>
#include <Colour.h>
#include <time.h>

#include <ModelTriangle.h>

#include <unordered_map>

#include<cstring>

#define WIDTH 320
#define HEIGHT 240



std::unordered_map<std::string,Colour> readColours(const std::string &filename){

	std::ifstream inputStream(filename);
	std::string line;
	std:: unordered_map<std::string,Colour> colourMap;

	while(inputStream && std::getline(inputStream,line)){

		std::string colourName = line;
		std::vector<std::string> splitColourName = split(colourName,' ');
		std::getline(inputStream,line);
		std::vector<std::string> splitLine = split(line,' ');
		colourMap[splitColourName[1]] =  Colour(stof(splitLine[1]), stof(splitLine[2]), stof(splitLine[3]));
		std::getline(inputStream,line);
	}

	return colourMap;

}

glm::vec2 getCanvasIntersectionPoint(glm::vec3 cameraPos, float focalLength, glm::vec3 vertex){

	float u = focalLength * (vertex.x/vertex.z) + (WIDTH/2);
	float v = focalLength * (vertex.y/vertex.z) + (HEIGHT/2);

	return glm::vec2(u,v);
	
}
void pointCloudRenderer(glm::vec3 coordinates, DrawingWindow &window,float focalLength, glm::vec3 cameraPos){

	
		glm::vec3 scaling = glm::vec3(240,240,240);
		glm::vec2 p1 = getCanvasIntersectionPoint(cameraPos,focalLength,coordinates*scaling);
		 window.setPixelColour(p1.x, p1.y, 0xFFFFFFFF);

	}

std::vector<ModelTriangle> readOBJFile(const std::string &filename,float scalingFactor,std::unordered_map<std::string,Colour> colourMap,float focalLength, glm::vec3 cameraPos, DrawingWindow &window ){


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
			pointCloudRenderer(glm::vec3(x,y,z),window,focalLength,cameraPos);
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
 void drawLine(CanvasPoint from, CanvasPoint to, DrawingWindow &window, Colour colour){

		float xDiff = to.x - from.x;
		float yDiff = to.y - from.y;
		float noOfSteps = std::fmax(abs(xDiff), abs(yDiff));
		float xStepSize = xDiff/noOfSteps;
		float yStepSize = yDiff/noOfSteps;

		float red = colour.red;
		float green = colour.green;
		float blue = colour.blue;
		uint32_t pixelColour = (255 << 24) + (int(red)<<16) + (int(green)<<8 ) + (int(blue));
		
		for( float i = 0.0; i<noOfSteps; i++){

			float x = from.x + (xStepSize * i);
			float y = from.y + (yStepSize * i);
			window.setPixelColour(round(x),round(y),pixelColour);
		}

 }
void drawTriangle(CanvasTriangle tr, DrawingWindow &window, Colour colour){

	drawLine(tr.v0(),tr.v1(),window,colour);
	drawLine(tr.v0(),tr.v2(),window,colour);
	drawLine(tr.v1(),tr.v2(),window,colour);

}

void drawFilledTriangle(CanvasTriangle tr, DrawingWindow &window, Colour colour){

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

	drawTriangle(tr,window,colour);
	drawLine(right,newPoint,window,colour);
	
	for(float Y = top.y; Y < y; Y++){

		float leftX = top.x - (top.x - x)*(Y-top.y)/(y- top.y);
		float rightX = top.x - (top.x - right.x)*(Y - top.y)/(right.y - top.y);

		drawLine(CanvasPoint(leftX,Y),CanvasPoint(rightX,Y),window,colour);
	}

	for(float Y = newPoint.y; Y < bottom.y; Y++){

		float leftX = newPoint.x - (newPoint.x - bottom.x)*(Y-newPoint.y)/(bottom.y-newPoint.y);
		float rightX = right.x - (right.x - bottom.x)*(Y - right.y)/(bottom.y - right.y);
		drawLine(CanvasPoint(leftX,Y),CanvasPoint(rightX,Y),window,colour);
	}


	drawTriangle(tr,window,Colour(255,255,255));	
}

void drawTexturedTriangle(CanvasTriangle canvasTriangle, CanvasTriangle textureTriangle,DrawingWindow &window, TextureMap textureMap){


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

		drawTriangle(canvasTriangle,window,Colour(255,255,255));	


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
			// if(i == WIDTH)
			// i = 0;
			 uint32_t colour = (255 << 24) + (int(rowInterpol.at(x).x) << 16) + ((int(rowInterpol.at(x).y)) << 8) + (int(rowInterpol.at(x).z));
			// i++;

			
			 window.setPixelColour(x, y, colour);


		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	} else if(event.key.keysym.sym == SDLK_u){

		CanvasPoint p1(rand() % WIDTH, rand() % HEIGHT);
		CanvasPoint p2(rand() % WIDTH, rand() % HEIGHT);
		CanvasPoint p3(rand() % WIDTH, rand() % HEIGHT);
		CanvasTriangle tr(p1,p2,p3);
		Colour colour(rand() % 256, rand() % 256, rand() % 256);
		drawTriangle(tr, window, colour);
	} else if(event.key.keysym.sym == SDLK_f){

		CanvasPoint p1(rand() % WIDTH, rand() % HEIGHT);
		CanvasPoint p2(rand() % WIDTH, rand() % HEIGHT);
		CanvasPoint p3(rand() % WIDTH, rand() % HEIGHT);
		CanvasTriangle tr(p1,p2,p3);
		Colour colour(rand() % 256, rand() % 256, rand() % 256);

		drawFilledTriangle(tr,window,colour);
	}
}



int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	srand(time(NULL));


	glm::vec3 cameraPos = glm::vec3(0.0,0.0,4.0);
	TextureMap map = TextureMap("src/texture.ppm");
	std::unordered_map<std::string,Colour> colourMap = readColours("src/cornell-box.mtl");
	std::vector<ModelTriangle> test = readOBJFile("src/cornell-box.obj",0.35,colourMap,2,cameraPos,window);


	// for(ModelTriangle mtr: test){
	// 	std::cout<<mtr<<" "<<mtr.colour<<" "<<std::endl;
	// }
	//  for (auto x : colourMap){

    // std::cout << x.first << " " <<
    //         x.second << std::endl;
	// 		}	
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		//draw(window);
		//pointCloudRenderer(test,window,2,cameraPos);
		// drawTexturedTriangle(CanvasTriangle(CanvasPoint(160, 10),CanvasPoint(300, 230),CanvasPoint(10, 150)),
		// CanvasTriangle(CanvasPoint(195, 5),CanvasPoint(395, 380),CanvasPoint(65, 330)),window,map);
		

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();

}

}

