#version 460 core
out vec4 FragColor;
//移除了Material的环境光颜色向量ambient，因为环境光颜色ambient在几乎所有情况下都等于漫反射颜色diffuse
struct Material {
    sampler2D diffuse;
    sampler2D specular;
    float shininess;
}; 

struct Light {
    // vec3 direction;
    vec3 position;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    float constant;
    float linear;
    float quadratic;
};

uniform vec3 viewPos;//观察者/摄像机的位置
uniform Material material;
uniform Light light;

in vec3 Normal;
in vec3 worldCoord;
in vec2 TexCoords;

void main()
{
    float distance    = length(light.position - worldCoord);
    float attenuation = 1.0 / (light.constant + light.linear * distance + 
                    light.quadratic * (distance * distance));

    //环境光照
    vec3 ambient = light.ambient * vec3(texture(material.diffuse, TexCoords));
    ambient *= attenuation;

    //漫反射光照
    vec3 norm = normalize(Normal);
    // vec3 lightDir = normalize(-light.direction);//平行光
    vec3 lightDir = normalize(light.position - worldCoord);//点光源
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = light.diffuse * diff * vec3(texture(material.diffuse, TexCoords));
    diffuse *= attenuation;

    //镜面光照
    vec3 viewDir = normalize(viewPos - worldCoord);//从片段指向摄像机
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);//32是高光的反光度(Shininess)，一个物体的反光度越高，反射光的能力越强，散射得越少，高光点就会越小。
    vec3 specular = spec * light.specular * vec3(texture(material.specular, TexCoords));
    specular *= attenuation;

    vec3 result = ambient + diffuse + specular;
    FragColor = vec4(result, 1.0);
}