#version 460 core
out vec4 FragColor;

void main()
{
    //不变的常量白色，保证了灯的颜色一直是亮的
    FragColor = vec4(1.0); // set alle 4 vector values to 1.0
}