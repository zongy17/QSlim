# 运行说明

Windows系统下使用CMake和MinGW-w64完成配置和编译

## 步骤
0. 确定安装有cmake和MingGW-w64，并将相应bin目录加入环境变量

<!-- 1. 阅读CMakeLists.txt修改头文件目录和库文件目录，将它们指向上传的压缩包内的openGL4.6下的对应路径 -->

2. 在项目目录下，依次执行下列命令完成配置和编译
```
mkdir build
cd build
cmake.exe -G "MinGW Makefiles" ..
mingw32-make.exe
```

3. 拷贝lib目录中的glfw.dll到build目录下（windows下使用动态链接库似乎必须要在运行目录下有对应动态库）

4. 运行.exe文件，注意命令行参数的选取。所有的.obj文件均在resource目录下，第一个参数为输入的.obj文件名，但输入时不要带".obj"后缀；第二个参数为网格简化的级数（simplification level）。resource目录下所有可用的.obj文件，即对应输入参数，如下列出：
- Arma.obj      => Arma
- Buddha.obj    => Buddha
- bunny.obj     => bunny
- cow.obj       => cow
- dinosaur.obj  => dinosaur
- dragon.obj    => dragon
- eight.uniform.obj => eight.uniform
- rocker.obj    => rocker
- sphere.obj    => sphere
- teddy.obj     => teddy
    
    类似FPS游戏，在程序中“W”“A”“S”“D”分别控制摄像机前、左、后、右平移。鼠标（光标）会被捕捉，移动鼠标会变换视角。鼠标滚轮可以缩放视域。

    在程序中“←”“→”方向键控制切换级数（level），长按“→”做进一步的简化，长按“←”返回上一级简化的网络。

```
B7.exe ../resource/cow 6 # 以cow.obj为输入，最多做6次简化
B7.exe ../resource/Arma 4 # 以Arma.obj为输入，最多做4次简化
```

