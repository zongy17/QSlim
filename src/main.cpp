#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "shader_class/shader.h"
#include "camera_class/camera.h"
#include "glad/gldebug.h"

#include <set>
#include <iostream>
#include "zyMesh.h"

#include <windows.h>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void processInput(GLFWwindow *window);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

// camera
Camera camera(glm::vec3(2.0f, 2.0f, 1.0f));
// 鼠标状态
bool CURSOR_DISABLED = true;//是否选择cursor追随模式（FPS相机）
bool MOUSE_LEFT_ON = false;
bool MOUSE_RIGHT_ON = false;
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

// lighting 全局vec3变量来表示光源在场景的世界空间坐标中的位置
glm::vec3 lightPos(3.0f, 3.0f, 0.0f);

struct HalfedgeHandle_to_cmp {
	size_t e_id;//用来找到边
	size_t vid0, vid1;
	float cost;//代价，用来维护红黑树中的偏序关系
	HalfedgeHandle_to_cmp(size_t e_id, size_t v0, size_t v1, float c): e_id(e_id), cost(c) {
		vid0 = std::min(v0, v1);
		vid1 = std::max(v0, v1);
	}
};
struct cmp {
	bool operator()(struct HalfedgeHandle_to_cmp const & a, \
					struct HalfedgeHandle_to_cmp const & b) const {
		if (a.cost == b.cost) return a.e_id < b.e_id;
		else return a.cost < b.cost;
	}
};

const float invertible_threshold = 1e-4;//判定为不可逆的阈值
std::set<HalfedgeHandle_to_cmp, cmp> min_heap;

// assumes that mesh has been setup
void initQEM(zyMesh & mesh) {
	min_heap.clear();
	// Step 1：利用周围各个面的p计算每个顶点的误差矩阵Q(面法向量和方程的计算在初始化载入时就可以做了)
	for (size_t vh = 0; vh < mesh.vertexList.size(); vh++) {
		mat44 Q = mat44(0.0f);
		//一个顶点的所有incomin或outgoing半边就对应了这个顶点周围的所有面
		std::vector<size_t> const & in_he_ids = mesh.vertexList[vh].incoming_halfedge_ids;
		for (unsigned vh_hei = 0; vh_hei < in_he_ids.size(); vh_hei++) {
			HalfedgeHandle const & heh = mesh.halfedgeList[in_he_ids[vh_hei]];
			// CVecf4 p = mesh.faces_p[heh.face_id];
			vec4 const & p = mesh.faceData[heh.face_id].Normal;
			Q += glm::outerProduct(p, p);
		}
		mesh.vertexData[vh].Q_vert = Q;
	}
	// Step 2：计算每条边进行缩而并产生新节点的Q_bar
	for (size_t e_id = 0; e_id < mesh.edgeList.size(); e_id++) {		
		HalfedgeHandle const & heh = mesh.halfedgeList[mesh.edgeList[e_id].halfedge_1_id];

		size_t v1_h = heh.from_vertex_id, v2_h = heh.to_vertex_id;
		mat44 Q1 = mesh.vertexData[v1_h].Q_vert, Q2 = mesh.vertexData[v2_h].Q_vert;
		mat44 Q_bar = Q1 + Q2;
		mat44 Q_bar_zero_LR = Q_bar;
		Q_bar_zero_LR[0][3] = 0.0f;// 最后一行（Last Row）赋值为0 0 0 1
		Q_bar_zero_LR[1][3] = 0.0f;
		Q_bar_zero_LR[2][3] = 0.0f;
		Q_bar_zero_LR[3][3] = 1.0f;
		vec4 v_bar;//合并后的顶点
		float cost;//合并的代价
		if (glm::determinant(Q_bar_zero_LR) > invertible_threshold) {// 可逆
			v_bar = glm::inverse(Q_bar_zero_LR) * vec4(0, 0, 0, 1);
			cost = glm::dot(v_bar, Q_bar * v_bar);
		} else {// 不可逆 => ﬁnd the optimal vertex along the segment v1v2
			// v_bar = k*v1 + (1-k)*v2
			// c = v_bar^T * Q_bar * v_bar = k^2*v1^T*Q_bar*v1 + (1-k)^2*v2^T*Q_bar*v2 + k*(1-k)*(v1^T*Q_bar*v2 + v2^T*Q_bar*v1)
			vec4 v1 = vec4(mesh.vertexData[v1_h].Position.x, mesh.vertexData[v1_h].Position.y, mesh.vertexData[v1_h].Position.z, 1);
			vec4 v2 = vec4(mesh.vertexData[v2_h].Position.x, mesh.vertexData[v2_h].Position.y, mesh.vertexData[v2_h].Position.z, 1);
			float v1Qv1 = glm::dot(v1, Q_bar * v1);
			float v1Qv2 = glm::dot(v1, Q_bar * v2);// 其实Q_bar是对称阵，因此v1Qv2==v2Qv1
			float v2Qv1 = glm::dot(v2, Q_bar * v1);
			float v2Qv2 = glm::dot(v2, Q_bar * v2);
			// dc/dk = 2*k*v1^T*Q_bar*v1 - 2*(1-k)*v2^T*Q_bar*v2  + (1-2*k)*(v1^T*Q_bar*v2 + v2^T*Q_bar*v1) = 0
			// 即 2k*c11 - 2*c22 + 2k*c22 + (c12+c21) - 2k*(c12+c21) = 0 即 k*c11 - c22 + k*c22 + 0.5*(c12+c21) - k*(c12+c21) = 0
			// 解出 k = (2*c22-c12-c21) / (2*c11+2*c22-2*c12-2*c21) = (c22 - 0.5*(c12+c21)) / (c11+c22-c12-c21)
			float k = (v2Qv2 - 0.5*(v1Qv2+v2Qv1)) / (v1Qv1+v2Qv2-v1Qv2-v2Qv1);
			if (0 <= k && k <= 1) {
				v_bar = k*v1 + (1-k)*v2;
				cost = glm::dot(v_bar, Q_bar * v_bar);
			} else {// If this also fails, we fall back on choosing v¯ from amongst the endpoints and the midpoint.
				vec4 vmid = vec4(0.5)*(v1 + v2);
				float cmid = glm::dot(vmid, Q_bar * vmid);
				if (cmid < v1Qv1 && cmid < v2Qv2){// c1即为v1Qv1 c2即为v2Qv2
					v_bar = vmid;
					cost = cmid;
				} else if (v1Qv1 < v2Qv2) {
					v_bar = v1;
					cost = v1Qv1;
				} else {
					v_bar = v2;
					cost = v2Qv2;
				}
			}
		}
		mesh.edgeData[heh.edge_id].Q_bar = Q_bar;
		mesh.edgeData[heh.edge_id].v_bar = v_bar;
		mesh.edgeData[heh.edge_id].cost  = cost;
		min_heap.insert(HalfedgeHandle_to_cmp(e_id, v1_h, v2_h, cost));
	}
}

size_t runQEM(zyMesh & mesh) {
	size_t old_n_vertices = mesh.vertexList.size();
	size_t plan_num_del = old_n_vertices / 2;

	if (old_n_vertices - plan_num_del < 4) {//简化完只剩不足4个顶点
		std::cout << "cannot simplify anymore!" << std::endl << std::endl;
		return 0;
	}
	//看起来每次简化都要先跑一次initQEM()，因为上次简化完之后存储在set中的句柄信息已经失效了！
	initQEM(mesh);
	
	size_t real_num_del = 0;
	for (int i = 0; i < plan_num_del; i++) {
		HalfedgeHandle_to_cmp top_halfedge = *min_heap.begin();//堆顶是代价最小的边 删除后 更新与这条边相连的边的代价
		min_heap.erase(min_heap.begin());//从堆中删除
		HalfedgeHandle & heh = mesh.halfedgeList[mesh.edgeList[top_halfedge.e_id].halfedge_1_id];

		mat44 Q_bar_del = mesh.edgeData[heh.edge_id].Q_bar;
		vec4 v_bar_del = mesh.edgeData[heh.edge_id].v_bar;

		VertexHandle & v_from = mesh.vertexList[heh.from_vertex_id];
		VertexHandle & v_to   = mesh.vertexList[heh.to_vertex_id];

		if (mesh.is_collapse_OK(heh)) {//合并不会导致拓扑改变
			real_num_del++;
			// 先要从set中删掉原来与v_from和v_to相连的边的记录
			// v_from
			for (unsigned ngb = 0; ngb < v_from.incoming_halfedge_ids.size(); ngb++){
				size_t he_id_del = v_from.incoming_halfedge_ids[ngb];
				size_t e_id_del  = mesh.halfedgeList[he_id_del].edge_id;
				min_heap.erase(HalfedgeHandle_to_cmp(e_id_del, mesh.halfedgeList[he_id_del].from_vertex_id,\
									mesh.halfedgeList[he_id_del].to_vertex_id, mesh.edgeData[e_id_del].cost));
			}
			for (unsigned ngb = 0; ngb < v_from.outgoing_halfedge_ids.size(); ngb++){
				size_t he_id_del = v_from.outgoing_halfedge_ids[ngb];
				size_t e_id_del  = mesh.halfedgeList[he_id_del].edge_id;
				min_heap.erase(HalfedgeHandle_to_cmp(e_id_del, mesh.halfedgeList[he_id_del].from_vertex_id,\
									mesh.halfedgeList[he_id_del].to_vertex_id, mesh.edgeData[e_id_del].cost));
			}
			// v_to
			for (unsigned ngb = 0; ngb < v_to.incoming_halfedge_ids.size(); ngb++){
				size_t he_id_del = v_to.incoming_halfedge_ids[ngb];
				size_t e_id_del  = mesh.halfedgeList[he_id_del].edge_id;
				min_heap.erase(HalfedgeHandle_to_cmp(e_id_del, mesh.halfedgeList[he_id_del].from_vertex_id,\
									mesh.halfedgeList[he_id_del].to_vertex_id, mesh.edgeData[e_id_del].cost));
			}
			for (unsigned ngb = 0; ngb < v_to.outgoing_halfedge_ids.size(); ngb++){
				size_t he_id_del = v_to.outgoing_halfedge_ids[ngb];
				size_t e_id_del  = mesh.halfedgeList[he_id_del].edge_id;
				min_heap.erase(HalfedgeHandle_to_cmp(e_id_del, mesh.halfedgeList[he_id_del].from_vertex_id,\
									mesh.halfedgeList[he_id_del].to_vertex_id, mesh.edgeData[e_id_del].cost));
			}
			
			// 将新点修改在v（collapse的终点）上
			mesh.zyCollapse(heh);
			mesh.vertexData[v_to.vertex_id].Position = vec3(v_bar_del.x, v_bar_del.y, v_bar_del.z);//修改终点的值为该边压缩之后的Vbar的坐标
			mesh.vertexData[v_to.vertex_id].Q_vert = Q_bar_del;			
			// 开始更新与最小代价边相连的边 也就是更新与新merged顶点相连的边
			for (unsigned i = 0; i < v_to.outgoing_halfedge_ids.size(); i++){
				HalfedgeHandle const & heh_v_to_outgoing = mesh.halfedgeList[v_to.outgoing_halfedge_ids[i]];
				size_t v_target = heh_v_to_outgoing.to_vertex_id;
				mat44 Q_bar_update = mesh.vertexData[v_target].Q_vert + mesh.vertexData[v_to.vertex_id].Q_vert;
				mat44 Q_bar_zero_LR = Q_bar_update;
				Q_bar_zero_LR[0][3] = 0.0f;// 最后一行（Last Row）赋值为0 0 0 1
				Q_bar_zero_LR[1][3] = 0.0f;
				Q_bar_zero_LR[2][3] = 0.0f;
				Q_bar_zero_LR[3][3] = 1.0f;

				vec4 v_bar_update;//合并后的顶点
				float cost_update;//合并的代价
				if (glm::determinant(Q_bar_zero_LR) > invertible_threshold) {// 可逆
					v_bar_update = glm::inverse(Q_bar_zero_LR) * vec4(0, 0, 0, 1);
					cost_update = glm::dot(v_bar_update, Q_bar_update * v_bar_update);
				} else {// 不可逆 => ﬁnd the optimal vertex along the segment v1v2
					// v_bar = k*v1 + (1-k)*v2
					// c = v_bar^T * Q_bar * v_bar = k^2*v1^T*Q_bar*v1 + (1-k)^2*v2^T*Q_bar*v2 + k*(1-k)*(v1^T*Q_bar*v2 + v2^T*Q_bar*v1)
					vec4 v1 = vec4(mesh.vertexData[v_target].Position.x, mesh.vertexData[v_target].Position.y, mesh.vertexData[v_target].Position.z, 1);
					vec4 v2 = vec4(mesh.vertexData[v_to.vertex_id].Position.x, mesh.vertexData[v_to.vertex_id].Position.y, mesh.vertexData[v_to.vertex_id].Position.z, 1);
					float v1Qv1 = glm::dot(v1, Q_bar_update * v1);
					float v1Qv2 = glm::dot(v1, Q_bar_update * v2);// 其实Q_bar是对称阵，因此v1Qv2==v2Qv1
					float v2Qv1 = glm::dot(v2, Q_bar_update * v1);
					float v2Qv2 = glm::dot(v2, Q_bar_update * v2);
					// dc/dk = 2*k*v1^T*Q_bar*v1 - 2*(1-k)*v2^T*Q_bar*v2  + (1-2*k)*(v1^T*Q_bar*v2 + v2^T*Q_bar*v1) = 0
					// 即 2k*c11 - 2*c22 + 2k*c22 + (c12+c21) - 2k*(c12+c21) = 0 即 k*c11 - c22 + k*c22 + 0.5*(c12+c21) - k*(c12+c21) = 0
					// 解出 k = (2*c22-c12-c21) / (2*c11+2*c22-2*c12-2*c21) = (c22 - 0.5*(c12+c21)) / (c11+c22-c12-c21)
					float k = (v2Qv2 - 0.5*(v1Qv2+v2Qv1)) / (v1Qv1+v2Qv2-v1Qv2-v2Qv1);
					if (0 <= k && k <= 1) {
						v_bar_update = k*v1 + (1-k)*v2;
						cost_update = glm::dot(v_bar_update, Q_bar_update * v_bar_update);
					} else {// If this also fails, we fall back on choosing v¯ from amongst the endpoints and the midpoint.
						vec4 vmid = vec4(0.5)*(v1 + v2);
						float cmid = glm::dot(vmid, Q_bar_update * vmid);
						if (cmid < v1Qv1 && cmid < v2Qv2){// c1即为v1Qv1 c2即为v2Qv2
							v_bar_update = vmid;
							cost_update = cmid;
						} else if (v1Qv1 < v2Qv2) {
							v_bar_update = v1;
							cost_update = v1Qv1;
						} else {
							v_bar_update = v2;
							cost_update = v2Qv2;
						}
					}
				}
				mesh.edgeData[heh_v_to_outgoing.edge_id].Q_bar = Q_bar_update;
				mesh.edgeData[heh_v_to_outgoing.edge_id].v_bar = v_bar_update;
				mesh.edgeData[heh_v_to_outgoing.edge_id].cost  = cost_update;
				assert(v_target != v_to.vertex_id);
				min_heap.insert(HalfedgeHandle_to_cmp(heh_v_to_outgoing.edge_id, v_target, v_to.vertex_id, cost_update));//重新加入新点与它周围顶点相连的边
			}
			for (unsigned i = 0; i < v_to.incoming_halfedge_ids.size(); i++){
				HalfedgeHandle const & heh_v_to_incoming = mesh.halfedgeList[v_to.incoming_halfedge_ids[i]];
				size_t v_src = heh_v_to_incoming.from_vertex_id;
				mat44 Q_bar_update = mesh.vertexData[v_src].Q_vert + mesh.vertexData[v_to.vertex_id].Q_vert;
				mat44 Q_bar_zero_LR = Q_bar_update;
				Q_bar_zero_LR[0][3] = 0.0f;// 最后一行（Last Row）赋值为0 0 0 1
				Q_bar_zero_LR[1][3] = 0.0f;
				Q_bar_zero_LR[2][3] = 0.0f;
				Q_bar_zero_LR[3][3] = 1.0f;

				vec4 v_bar_update;//合并后的顶点
				float cost_update;//合并的代价
				if (glm::determinant(Q_bar_zero_LR) > invertible_threshold) {// 可逆
					v_bar_update = glm::inverse(Q_bar_zero_LR) * vec4(0, 0, 0, 1);
					cost_update = glm::dot(v_bar_update, Q_bar_update * v_bar_update);
				} else {// 不可逆 => ﬁnd the optimal vertex along the segment v1v2
					// v_bar = k*v1 + (1-k)*v2
					// c = v_bar^T * Q_bar * v_bar = k^2*v1^T*Q_bar*v1 + (1-k)^2*v2^T*Q_bar*v2 + k*(1-k)*(v1^T*Q_bar*v2 + v2^T*Q_bar*v1)
					vec4 v1 = vec4(mesh.vertexData[v_src].Position.x, mesh.vertexData[v_src].Position.y, mesh.vertexData[v_src].Position.z, 1);
					vec4 v2 = vec4(mesh.vertexData[v_to.vertex_id].Position.x, mesh.vertexData[v_to.vertex_id].Position.y, mesh.vertexData[v_to.vertex_id].Position.z, 1);
					float v1Qv1 = glm::dot(v1, Q_bar_update * v1);
					float v1Qv2 = glm::dot(v1, Q_bar_update * v2);// 其实Q_bar是对称阵，因此v1Qv2==v2Qv1
					float v2Qv1 = glm::dot(v2, Q_bar_update * v1);
					float v2Qv2 = glm::dot(v2, Q_bar_update * v2);
					// dc/dk = 2*k*v1^T*Q_bar*v1 - 2*(1-k)*v2^T*Q_bar*v2  + (1-2*k)*(v1^T*Q_bar*v2 + v2^T*Q_bar*v1) = 0
					// 即 2k*c11 - 2*c22 + 2k*c22 + (c12+c21) - 2k*(c12+c21) = 0 即 k*c11 - c22 + k*c22 + 0.5*(c12+c21) - k*(c12+c21) = 0
					// 解出 k = (2*c22-c12-c21) / (2*c11+2*c22-2*c12-2*c21) = (c22 - 0.5*(c12+c21)) / (c11+c22-c12-c21)
					float k = (v2Qv2 - 0.5*(v1Qv2+v2Qv1)) / (v1Qv1+v2Qv2-v1Qv2-v2Qv1);
					if (0 <= k && k <= 1) {
						v_bar_update = k*v1 + (1-k)*v2;
						cost_update = glm::dot(v_bar_update, Q_bar_update * v_bar_update);
					} else {// If this also fails, we fall back on choosing v¯ from amongst the endpoints and the midpoint.
						vec4 vmid = vec4(0.5)*(v1 + v2);
						float cmid = glm::dot(vmid, Q_bar_update * vmid);
						if (cmid < v1Qv1 && cmid < v2Qv2){// c1即为v1Qv1 c2即为v2Qv2
							v_bar_update = vmid;
							cost_update = cmid;
						} else if (v1Qv1 < v2Qv2) {
							v_bar_update = v1;
							cost_update = v1Qv1;
						} else {
							v_bar_update = v2;
							cost_update = v2Qv2;
						}
					}
				}
				mesh.edgeData[heh_v_to_incoming.edge_id].Q_bar = Q_bar_update;
				mesh.edgeData[heh_v_to_incoming.edge_id].v_bar = v_bar_update;
				mesh.edgeData[heh_v_to_incoming.edge_id].cost  = cost_update;
				assert(v_src != v_to.vertex_id);
				min_heap.insert(HalfedgeHandle_to_cmp(heh_v_to_incoming.edge_id, v_src, v_to.vertex_id, cost_update));//重新加入新点与它周围顶点相连的边
			}
		}
	}

	return real_num_del;
}

GLFWwindow * initGLFW(unsigned int const major_version, unsigned int const minor_version) {
    // glfw: initialize and configure
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);//如果需要请求调试输出

    // glfw window creation
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "B-7", NULL, NULL);
    if (window == NULL){
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    if (CURSOR_DISABLED) {
        glfwSetCursorPosCallback(window, mouse_callback);
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);// tell GLFW to capture our mouse
    }
    glfwSetScrollCallback(window, scroll_callback);
    return window;
}

void doSimplify(zyMesh & mesh, std::string filename_nopfx, std::vector<zyMesh> & meshList) {
	printf("enter doSimplify...\n");
	unsigned long start_time, end_time;
    start_time = GetTickCount();
    size_t numVerts_del = runQEM(mesh);// 做完简化之后的mesh已经不是原来的样子了（vertexData已变）
    end_time = GetTickCount();
    std::cout << "level " << meshList.size() << " simplification cost: " << end_time-start_time << " ms" \
			  << " with " << numVerts_del << " vertices deleted."  << std::endl;
    std::string filename = filename_nopfx + "-smpf-" + std::to_string(meshList.size()) + ".obj";
	std::cout << "write to " << filename << std::endl;
	mesh.write_obj(filename.c_str());// 先写出，再读入
	// std::cout << "load mesh: " << filename << std::endl;
	meshList.push_back(zyMesh(filename, true, true));//只有新加入的这一个mesh是可以进一步操作的
	printf("leave doSimplify...\n");
}

void genMeshes(std::string filename_nopfx, unsigned level, std::vector<zyMesh> & meshList) {
	meshList.push_back(zyMesh(filename_nopfx+".obj", true, true));
	zyMesh tmp_mesh;
	unsigned l = 0; unsigned long start_time, end_time;
	std::string in_filename = filename_nopfx+".obj", out_filename;
	while (l < level) {
		l++;
		out_filename = filename_nopfx + "-smpf-" + std::to_string(l) + ".obj";
		
		printf("doing %d-th level of simplification...  ", l);
		start_time = GetTickCount();
		tmp_mesh = zyMesh(in_filename, true, false);//中间网格不需要设置openGL相关的数据
		size_t numVerts_del = runQEM(tmp_mesh);// 做完简化之后的mesh已经不是原来的样子了（vertexData已变）
		tmp_mesh.write_obj(out_filename.c_str());
		meshList.push_back(zyMesh(out_filename, true, true));
		end_time = GetTickCount();
		printf("cost %u ms with %u vertices deleted.\n", end_time-start_time, numVerts_del);

		in_filename = out_filename;
	}
}

int smpf_level = 0, max_level = 3;
float smpf_level_f = (float)smpf_level;
void renderMain(GLFWwindow * window, std::string filename_nopfx, Shader & ourShader, std::vector<zyMesh> & meshList) {
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        processInput(window);

		// if (meshList.size()-1 < max_level && meshList.size()-1 <= smpf_level_f) {// 如果当前meshList内已经有足够的（max_level）网格，则不再细分
		// 	zyMesh tmp_mesh(meshList[smpf_level]);
        //     doSimplify(tmp_mesh, filename_nopfx, meshList);
        // }

        // render
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // don't forget to enable shader before setting uniforms
        ourShader.use();

        // view/projection transformations
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();
        ourShader.setMat4("projection", projection);
        ourShader.setMat4("view", view);

        // render the loaded model
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, glm::vec3(0.0f, 0.0f, 0.0f)); // translate it down so it's at the center of the scene
		if (filename_nopfx.find("bunny") == std::string::npos)
			model = glm::scale(model, glm::vec3(0.02f, 0.02f, 0.02f));	// it's a bit too big for our scene, so scale it down
		else
			model = glm::scale(model, glm::vec3(30.0f, 30.0f, 30.0f));
        ourShader.setMat4("model", model);

        ourShader.setBool("sgColor", true);
        ourShader.setVec3("backColor", glm::vec3(160.0f, 32.0f, 240.0f)/255.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//面
        // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//边（网格线 wireframe）
        // glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);//点
		// int ptr = (subdivision_level > meshList.size()-1) ? meshList.size()-1 : subdivision_level;
        meshList[smpf_level].Draw(ourShader);

        // 面+边（点）模式，则先画面，再画边（点）
        ourShader.setBool("sgColor", true);
        ourShader.setVec3("backColor", glm::vec3(0.0f, 0.0f, 0.0f));
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);//在面的基础上叠加上边
        glEnable(GL_POLYGON_OFFSET_LINE);
        // glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);//在面的基础上叠加上点
        // glEnable(GL_POLYGON_OFFSET_POINT);
        glPolygonOffset(-0.5, -0.5); // move closer to camera
        meshList[smpf_level].Draw(ourShader);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

int main(int argc, char* argv[])
{
    GLFWwindow* window = initGLFW(4, 6);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){// glad: load all OpenGL function pointers
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    glEnable(GL_DEPTH_TEST);// configure global opengl state
    glLineWidth(1.0f);
    glPointSize(2.0f);

    glEnable(GL_CULL_FACE);//面剔除，如果想要背面也画出来（尤其是GL_LINE和GL_POINT模式下），则关掉这两行
    glCullFace(GL_BACK);//剔除掉背向面
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    // build and compile shaders
    Shader ourShader("../resource/model_loading.vs", "../resource/model_loading.fs");

    std::string filename_nopfx = std::string(argv[1]);
	max_level = atoi(argv[2]);
    
	std::vector<zyMesh> meshList;
	genMeshes(filename_nopfx, max_level, meshList);
	// meshList.push_back(zyMesh(filename_nopfx+".obj", true, true));

    // render loop
    renderMain(window, filename_nopfx, ourShader, meshList);
    
    // glfw: terminate, clearing all previously allocated GLFW resources.
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);


        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {//正在被按下
            if (MOUSE_LEFT_ON == false) {//之前没被按下
                MOUSE_LEFT_ON = true;
                //保存下来lastX和lastY
                double xpos, ypos;
                glfwGetCursorPos(window, &xpos, &ypos);
                lastX = xpos;
                lastY = ypos;
            } else {//之前一直被按着
                //重新获得当前的位置
                double xpos, ypos;
                glfwGetCursorPos(window, &xpos, &ypos);
                float xoffset = xpos - lastX;
                float yoffset = lastY - ypos;
                lastX = xpos;
                lastY = ypos;

                camera.ProcessMouseMovement(-xoffset, -yoffset);//更改摄像机位置
            }
        } else {//被松开
            MOUSE_LEFT_ON = false;
        }

        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {//正在被按下
            if (MOUSE_RIGHT_ON == false) {//之前没被按下
                MOUSE_RIGHT_ON = true;
                //保存下来lastX和lastY
                double xpos, ypos;
                glfwGetCursorPos(window, &xpos, &ypos);
                lastX = xpos;
                lastY = ypos;
            } else {//之前一直被按着
                //重新获得当前的位置
                double xpos, ypos;
                glfwGetCursorPos(window, &xpos, &ypos);
                float sensitivity = 3e-3;
                float xoffset = (xpos - lastX)*sensitivity;
                float yoffset = (lastY - ypos)*sensitivity;
                lastX = xpos;
                lastY = ypos;

                camera.ProcessKeyboard(RIGHT, -xoffset);
                camera.ProcessKeyboard(UP   , -yoffset);
            }
        } else {//被松开
            MOUSE_RIGHT_ON = false;
        }
    
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS){
        smpf_level_f = std::max(0.0f, smpf_level_f-0.01f);
        smpf_level = (int)smpf_level_f;
    }
        
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS){
        smpf_level_f = std::min((float)max_level, smpf_level_f+0.01f);
        smpf_level = (int)smpf_level_f;
    }
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}