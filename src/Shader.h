#ifndef SHADER_H
#define SHADER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <GL/glew.h>

class Shader {
public:
	GLuint Program;

	Shader(const GLchar* vertexPath, const GLchar* fragmentPath)
	{
		std::string vertexCode;
		std::string fragmentCode;

		std::ifstream vShaderFile;
		std::ifstream fShaderFile;

		vShaderFile.exceptions (std::ifstream::badbit);
		fShaderFile.exceptions (std::ifstream::badbit);

		try
		{
			vShaderFile.open(vertexPath);
			fShaderFile.open(fragmentPath);

			std::stringstream vShaderStream, gShaderStream, fShaderStream;

			vShaderStream << vShaderFile.rdbuf();
			fShaderStream << fShaderFile.rdbuf();

			vShaderFile.close();
			fShaderFile.close();

			vertexCode = vShaderStream.str();
			fragmentCode = fShaderStream.str();
		}
		catch (std::ifstream::failure e) { std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl; }

		const GLchar * vShaderCode = vertexCode.c_str();
		const GLchar * fShaderCode = fragmentCode.c_str();

		GLuint vertex, geometry, fragment;
		GLint success;
		GLchar infoLog[512];
		
		
	
		vertex = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vertex, 1, &vShaderCode, NULL);
		glCompileShader(vertex);
		glGetShaderiv(vertex, GL_COMPILE_STATUS, &success);
				
		
		if (!success)
		{
			std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
		}

		fragment = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fragment, 1, &fShaderCode, NULL);
		glCompileShader(fragment);
		glGetShaderiv(fragment, GL_COMPILE_STATUS, &success);
		if (!success)
		{
			std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
		}


		this->Program = glCreateProgram();
		glAttachShader(this->Program, vertex);
		glAttachShader(this->Program, fragment);
		glLinkProgram(this->Program);

		glDeleteShader(vertex);
		glDeleteShader(fragment);
	}


	void Use()
	{
		glUseProgram(this->Program);
	}
};

#endif
