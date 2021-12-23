#version 460 core
out vec4 FragColor;
  
uniform vec3 objectColor;
uniform vec3 lightColor;
uniform vec3 lightPos;//光源的位置是一个静态变量
uniform vec3 viewPos;//观察者/摄像机的位置

in vec3 Normal;
in vec3 worldCoord;

void main()
{
    //环境光照
    float ambientStrength = 0.1;
    vec3 ambient = ambientStrength * lightColor;
    //漫反射光照
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - worldCoord);//从片段指向光源
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;
    //镜面光照
    float specularStrength = 1.0;
    vec3 viewDir = normalize(viewPos - worldCoord);//从片段指向摄像机
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);//32是高光的反光度(Shininess)
    //一个物体的反光度越高，反射光的能力越强，散射得越少，高光点就会越小。
    vec3 specular = specularStrength * spec * lightColor;

    vec3 result = (ambient + diffuse + specular) * objectColor;
    FragColor = vec4(result, 1.0);
    // FragColor = vec4(lightColor * objectColor, 1.0);
}