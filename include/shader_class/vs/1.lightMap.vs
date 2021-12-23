#version 460 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 Normal;
out vec3 worldCoord;//在世界空间中进行所有的光照计算，因此我们需要一个在世界空间中的顶点位置
out vec2 TexCoords;

void main()
{
    worldCoord = vec3(model * vec4(aPos, 1.0));
    
    //Normal = aNormal;
    //法线矩阵：模型矩阵左上角的逆矩阵的转置矩阵，用它来移除对法向量错误缩放的影响
    //Normal = mat3(transpose(inverse(model))) * aNormal;
    Normal = transpose(inverse(mat3(model))) * aNormal;

    gl_Position = projection * view * model * vec4(aPos, 1.0);

    TexCoords = aTexCoords;
}