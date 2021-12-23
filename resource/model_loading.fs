#version 460 core
out vec4 FragColor;

in vec3 Normal;
in vec2 TexCoords;

uniform sampler2D texture_diffuse1;
uniform bool sgColor;
uniform vec3 backColor;

void main()
{    
    if (sgColor == false)//正常纹理
        // FragColor = texture(texture_diffuse1, TexCoords);
        FragColor = vec4(Normal, 1.0f);
    else //单色
        FragColor = vec4(backColor, 1.0f);

}