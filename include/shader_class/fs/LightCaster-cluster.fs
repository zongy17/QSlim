#version 460 core
out vec4 FragColor;
//移除了Material的环境光颜色向量ambient，因为环境光颜色ambient在几乎所有情况下都等于漫反射颜色diffuse
struct Material {
    sampler2D diffuse;
    sampler2D specular;
    float shininess;
}; 

struct Light {
    vec3 direction;//聚光同时需要位置和方向向量
    vec3 position;
    float cutOff;//切光角
    float outterCutOff;

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
    // float distance    = length(light.position - worldCoord);
    // float attenuation = 1.0 / (light.constant + light.linear * distance + 
    //                 light.quadratic * (distance * distance));
    vec3 lightDir = normalize(light.position - worldCoord);
    float theta = dot(lightDir, normalize(-light.direction));
    float epsilon = light.cutOff - light.outterCutOff;
    float intensity = clamp((theta - light.outterCutOff) / epsilon, 0.0, 1.0);    

        // ambient
        vec3 ambient = light.ambient * texture(material.diffuse, TexCoords).rgb;
        
        // diffuse 
        vec3 norm = normalize(Normal);
        float diff = max(dot(norm, lightDir), 0.0);
        vec3 diffuse = light.diffuse * diff * texture(material.diffuse, TexCoords).rgb;  
        
        // specular
        vec3 viewDir = normalize(viewPos - worldCoord);
        vec3 reflectDir = reflect(-lightDir, norm);  
        float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);
        vec3 specular = light.specular * spec * texture(material.specular, TexCoords).rgb;  
        
        // (soft edges)
        diffuse  *= intensity;
        specular *= intensity;

        // attenuation
        float distance    = length(light.position - worldCoord);
        float attenuation = 1.0 / (light.constant + light.linear * distance + light.quadratic * (distance * distance));    
        ambient  *= attenuation; 
        diffuse   *= attenuation;
        specular *= attenuation;   
            
        vec3 result = ambient + diffuse + specular;
        FragColor = vec4(result, 1.0);
}