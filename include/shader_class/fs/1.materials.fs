#version 460 core
out vec4 FragColor;
  
struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;    
    float shininess;
}; 

struct Light {
    vec3 position;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

uniform vec3 viewPos;//观察者/摄像机的位置
uniform Material material;
uniform Light light;

in vec3 Normal;
in vec3 worldCoord;

void main()
{
    //环境光照
    vec3 ambient = material.ambient * light.ambient;

    //漫反射光照
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(light.position - worldCoord);//从片段指向光源
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = light.diffuse * (diff * material.diffuse);

    //镜面光照
    vec3 viewDir = normalize(viewPos - worldCoord);//从片段指向摄像机
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);//32是高光的反光度(Shininess)，一个物体的反光度越高，反射光的能力越强，散射得越少，高光点就会越小。
    vec3 specular = spec * light.specular * material.specular;

    vec3 result = ambient + diffuse + specular;
    FragColor = vec4(result, 1.0);
}