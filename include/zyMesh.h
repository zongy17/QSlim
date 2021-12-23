#ifndef __ZY__MESH__
#define __ZY__MESH__

#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "shader_class/shader.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include <limits.h>
#include <unordered_set>


typedef glm::mat4x4 mat44;
typedef glm::vec3 vec3;
typedef glm::vec4 vec4;
typedef glm::vec2 vec2;

struct VertexData {
    // position
    vec3 Position;
    // normal
    vec3 Normal;
    // texCoords
    vec2 TexCoords;
    // for qslim
    mat44 Q_vert;
    VertexData() {Position = Normal = vec3(0.0f); TexCoords = vec2(0.0f); Q_vert = mat44(0.0f);}
};

struct FaceData {
    vec4 Normal;
    FaceData() {Normal = vec4(0.0f);}
};

struct EdgeData {
    //for qslim
    float cost;
    mat44 Q_bar;
    vec4 v_bar;
    EdgeData() {cost = 0.0; Q_bar = mat44(0.0f); v_bar = vec4(0.0f);}
};

#define IndexNotValid UINT64_MAX
class VertexHandle {
public:
    size_t vertex_id;
    size_t smpf_vid = IndexNotValid;//smpf_vid是经过删除部分网格后重整后的结点序号
    std::vector<size_t> outgoing_halfedge_ids;
    std::vector<size_t> incoming_halfedge_ids;
    bool deleted = false;
    VertexHandle() {outgoing_halfedge_ids.clear(); incoming_halfedge_ids.clear();}
    VertexHandle(size_t vid) { vertex_id = vid; outgoing_halfedge_ids.clear(); incoming_halfedge_ids.clear();}
};

class HalfedgeHandle {
public:
    // Halfedge* prev, * next, * opposive;
    size_t halfedge_id;
    size_t edge_id;
    size_t prev_halfedge_id, next_halfedge_id, oppo_halfedge_id;
    size_t to_vertex_id, from_vertex_id;
    size_t face_id;
    bool deleted = false;
    HalfedgeHandle() {halfedge_id = edge_id = prev_halfedge_id = next_halfedge_id = oppo_halfedge_id = to_vertex_id = from_vertex_id = face_id = IndexNotValid;}

    bool operator == (HalfedgeHandle const & b) const {
        return halfedge_id==b.halfedge_id && edge_id==b.edge_id && prev_halfedge_id==b.prev_halfedge_id && next_halfedge_id==b.next_halfedge_id &&\
               oppo_halfedge_id==b.oppo_halfedge_id && to_vertex_id==b.to_vertex_id && from_vertex_id==b.from_vertex_id && face_id==b.face_id;
    }
    bool operator != (HalfedgeHandle const & b) const {return !this->operator==(b); }
};

class EdgeHandle {
public:
    size_t halfedge_1_id, halfedge_2_id;
    size_t h1_to_vid, h1_from_vid;
    bool deleted = false;
    EdgeHandle() {halfedge_1_id = halfedge_2_id = IndexNotValid;}
    EdgeHandle(size_t h1_id): halfedge_1_id(h1_id), halfedge_2_id(IndexNotValid) {}
};

class FaceHandle {
public:
    size_t start_halfegde_id;
    std::vector<size_t> vertex_ids;
    bool deleted = false;
    FaceHandle() {start_halfegde_id = IndexNotValid; vertex_ids.clear();}
    FaceHandle(size_t s_he_id, std::vector<size_t> const & vids): start_halfegde_id(s_he_id), vertex_ids(vids) {}
};

struct Texture {
    unsigned int id;
    std::string type;
    std::string path;  // 我们储存纹理的路径用于与其它纹理进行比较
};

class zyMesh {
public:
    std::vector<VertexData> vertexData;
    std::vector<VertexHandle> vertexList;
    std::vector<EdgeData> edgeData;
    std::vector<EdgeHandle> edgeList;
    std::vector<FaceData> faceData;
    std::vector<FaceHandle> faceList;
    std::vector<HalfedgeHandle> halfedgeList;
    std::vector<unsigned int> indices;

    // used for OpenGL drawing
    unsigned int VAO;
    unsigned int VBO, EBO;
    std::vector<Texture>      textures;

    zyMesh() {
        vertexData.clear();
        vertexList.clear();
        edgeData.clear();
        edgeList.clear();
        faceData.clear();
        faceList.clear(); 
        halfedgeList.clear();
        indices.clear(); 
        textures.clear();
    }
    zyMesh(std::string filename, bool _need_to_calc_normal=false, bool _setup_mesh=true) {
        vertexData.clear();
        vertexList.clear();
        edgeData.clear();
        edgeList.clear();
        faceData.clear();
        faceList.clear(); 
        halfedgeList.clear();
        indices.clear(); 

        textures.clear();
        loadData(filename);
        if (_need_to_calc_normal)
            calcNormal();
        
        if (_setup_mesh)
            setupMesh();
    }
    zyMesh(const zyMesh & b) {
        textures = b.textures;
        VAO = b.VAO;
        VBO = b.VBO;
        EBO = b.EBO;
        vertexData = b.vertexData;
        vertexList = b.vertexList;
        faceData = b.faceData;
        faceList = b.faceList;
        indices = b.indices;
        edgeData = b.edgeData; 
        edgeList = b.edgeList;
        halfedgeList = b.halfedgeList;
    }

    void setupMesh() {
        // create buffers/arrays
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);
        // load data into vertex buffers
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        // A great thing about structs is that their memory layout is sequential for all its items.
        // The effect is that we can simply pass a pointer to the struct and it translates perfectly to a glm::vec3/2 array which
        // again translates to 3/2 floats which translates to a byte array.
        glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(VertexData), &vertexData[0], GL_STATIC_DRAW);  

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);

        // set the vertex attribute pointers
        // vertex Positions
        glEnableVertexAttribArray(0);	
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexData), (void*)0);
        // vertex normals
        glEnableVertexAttribArray(1);	
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(VertexData), (void*)offsetof(VertexData, Normal));
        // vertex texture coords
        glEnableVertexAttribArray(2);	
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(VertexData), (void*)offsetof(VertexData, TexCoords));

        glBindVertexArray(0);
    }

    // 计算面和顶点的法向量
    void calcNormal() {
        // 假定是三角形
        for (size_t f = 0; f < faceList.size(); f++) {
            size_t vids[3];
            size_t start_halfedge_id = faceList[f].start_halfegde_id;
            size_t halfedge_id = start_halfedge_id;
            unsigned i = 0;
            do {
                vids[i++] = halfedgeList[halfedge_id].from_vertex_id;
                halfedge_id = halfedgeList[halfedge_id].next_halfedge_id;
            } while (halfedge_id != start_halfedge_id);
            vec3 v0v1 = vertexData[vids[1]].Position - vertexData[vids[0]].Position;
            vec3 v0v2 = vertexData[vids[2]].Position - vertexData[vids[0]].Position;
            vec3 Vabc = glm::normalize(glm::cross(v0v1, v0v2)); //v0v1.cross(v0v2).stableNormalized();
            float a = Vabc.x, b = Vabc.y, c = Vabc.z;
            float d = - ( a*vertexData[vids[0]].Position.x \
                        + b*vertexData[vids[0]].Position.y \
                        + c*vertexData[vids[0]].Position.z);
            faceData[f].Normal = vec4(a, b, c, d);
        }
        for (size_t vh = 0; vh < vertexList.size(); vh++) {
            vec3 norm_wArea = vec3(0.0f);
            //一个顶点的所有incomin或outgoing半边就对应了这个顶点周围的所有面
            for (unsigned vh_hei = 0; vh_hei < vertexList[vh].incoming_halfedge_ids.size(); vh_hei++) {
                HalfedgeHandle const & heh = halfedgeList[vertexList[vh].incoming_halfedge_ids[vh_hei]];
                size_t face_id = heh.face_id;
                size_t v1_id = heh.from_vertex_id, v2_id = halfedgeList[heh.next_halfedge_id].to_vertex_id;
                vec3 v0v1 = vertexData[v1_id].Position - vertexData[vh].Position;
                vec3 v0v2 = vertexData[v2_id].Position - vertexData[vh].Position;
                float l1 = sqrt(v0v1.x*v0v1.x + v0v1.y*v0v1.y + v0v1.z*v0v1.z);
                float l2 = sqrt(v0v2.x*v0v2.x + v0v2.y*v0v2.y + v0v2.z*v0v2.z);
                float angle = acos(glm::dot(v0v1,v0v2) / (l1*l2));
                vec4 const & faceNormal = faceData[face_id].Normal;
                norm_wArea += vec3(faceNormal.x, faceNormal.y, faceNormal.z) * angle;
                // norm_wArea += faceNormals_wArea[face_id] * angle;
            }
            vertexData[vh].Normal = glm::normalize(norm_wArea);// norm_wArea.normalized();
        }
    }

    size_t add_vertex(vec3 const & p) {
        VertexData vd;
        vd.Position = p;
        size_t vid = vertexList.size();
        vertexData.push_back(vd);
        vertexList.push_back(VertexHandle(vid));
        return vid;
    }
    
    // make sure vertex_ids already exist in the mesh, and sequence in vertex_ids is anti-clockwise and point-outward
    size_t add_face(std::vector<size_t> const & vertex_ids) {
        size_t num_edges = vertex_ids.size();
        if (num_edges < 3) {
            std::cout << "Error face with less than 3 edges!" << std::endl;
            exit(1);
        }
        size_t face_id = faceList.size();
        size_t start_he_id =  halfedgeList.size();
        //有多少个顶点，就有多少条边
        for (unsigned ie = 0; ie < num_edges; ie++) {
            HalfedgeHandle heh;
            size_t self_id = halfedgeList.size();
            heh.halfedge_id = self_id;
            heh.from_vertex_id = vertex_ids[ie];
            heh.to_vertex_id   = vertex_ids[(ie+1)%num_edges];
            //根据from_vertex来找oppo_half_edge_id
            for (unsigned fr_vt_he_i = 0; \
                fr_vt_he_i < vertexList[heh.from_vertex_id].incoming_halfedge_ids.size(); fr_vt_he_i++) {
                size_t oppo_heh_id = vertexList[heh.from_vertex_id].incoming_halfedge_ids[fr_vt_he_i];
                HalfedgeHandle const & oppo_heh = halfedgeList[oppo_heh_id];
                if (oppo_heh.from_vertex_id == heh.to_vertex_id){
                    heh.oppo_halfedge_id = oppo_heh_id;
                    halfedgeList[oppo_heh_id].oppo_halfedge_id = self_id;
                    break;
                }
            }
            vertexList[heh.from_vertex_id].outgoing_halfedge_ids.push_back(self_id);//由点索引半边
            vertexList[heh.to_vertex_id].incoming_halfedge_ids.push_back(self_id);
            heh.face_id = face_id;
            heh.prev_halfedge_id = (ie==0) ? (self_id-1+num_edges) : (self_id-1);
            heh.next_halfedge_id = (ie==(num_edges-1)) ? (self_id+1-num_edges) : (self_id+1);

            if (heh.oppo_halfedge_id != IndexNotValid){// 确定了半边对应边的关系，对边已经把edge加过了
                heh.edge_id = halfedgeList[heh.oppo_halfedge_id].edge_id;
                edgeList[heh.edge_id].halfedge_2_id = self_id;
            } else {
                heh.edge_id = edgeList.size();
                EdgeHandle newegde(self_id);
                newegde.h1_from_vid = heh.from_vertex_id;
                newegde.h1_to_vid = heh.to_vertex_id;
                edgeList.push_back(newegde);
                edgeData.push_back(EdgeData());//先压一个空的
            }

            halfedgeList.push_back(heh);
            indices.push_back(vertex_ids[ie]);
        }
        faceList.push_back(FaceHandle(start_he_id, vertex_ids));
        faceData.push_back(FaceData());//先压一个空的
        
        return start_he_id;
    }

    float calcFaceArea(const HalfedgeHandle & _heh) {
        size_t vids[3];
        size_t start_halfedge_id = faceList[_heh.face_id].start_halfegde_id;
        size_t halfedge_id = start_halfedge_id;
        unsigned i = 0;
        do {
            vids[i++] = halfedgeList[halfedge_id].from_vertex_id;
            halfedge_id = halfedgeList[halfedge_id].next_halfedge_id;
        } while (halfedge_id != start_halfedge_id);
        vec3 v0v1 = vertexData[vids[1]].Position - vertexData[vids[0]].Position;
        vec3 v0v2 = vertexData[vids[2]].Position - vertexData[vids[0]].Position;
        vec3 crs_prd = glm::cross(v0v1, v0v2);
        return 0.5*sqrt(crs_prd.x*crs_prd.x + crs_prd.y*crs_prd.y + crs_prd.z*crs_prd.z);
    }

    std::vector<size_t> write_face(HalfedgeHandle & heh) {
        std::vector<size_t> vids;
        size_t start_he_id = heh.halfedge_id, he_id = start_he_id;
        do {
            HalfedgeHandle & _hh = halfedgeList[he_id];
            VertexHandle const & _vh = vertexList[_hh.from_vertex_id];
            vids.push_back(_vh.smpf_vid);
            he_id = _hh.next_halfedge_id;
            _hh.deleted = true;
        } while (he_id != start_he_id);
        return vids;
    }

    void write_obj(const char *  filename) {
        // default: only vertices and faces
        FILE * fid = fopen(filename, "w+");
        size_t output_id = 0;
        for (size_t v = 0; v < vertexList.size(); v++) {
            if (!vertexList[v].deleted){
                vertexList[v].smpf_vid = output_id++;
                // fprintf(fid, "v %f %f %f\n", vertexData[v].Position.x(), vertexData[v].Position.y(), vertexData[v].Position.z());
                fprintf(fid, "v %f %f %f\n", vertexData[v].Position.x, vertexData[v].Position.y, vertexData[v].Position.z);
            }
        }
        // 用半边来找面
        for (size_t he_id = 0; he_id < halfedgeList.size(); he_id++){
            if (!halfedgeList[he_id].deleted){
                std::vector<size_t> vids = write_face(halfedgeList[he_id]);
                fprintf(fid, "f");
                for (unsigned i = 0; i < vids.size(); i++)
                    fprintf(fid, " %u", vids[i]+1);//注意obj格式+1
                fprintf(fid, "\n");
            }
        }
        fclose(fid);
    }

    void Draw(Shader &shader) {
        // bind appropriate textures
        unsigned int diffuseNr  = 1;
        unsigned int specularNr = 1;
        unsigned int normalNr   = 1;
        unsigned int heightNr   = 1;
        unsigned int reflectionNr = 1;
        for(unsigned int i = 0; i < textures.size(); i++)
        {
            glActiveTexture(GL_TEXTURE0 + i); // active proper texture unit before binding
            // retrieve texture number (the N in diffuse_textureN)
            std::string number;
            std::string name = textures[i].type;
            if(name == "texture_diffuse")
                number = std::to_string(diffuseNr++);
            else if(name == "texture_specular")
                number = std::to_string(specularNr++); // transfer unsigned int to stream
            else if(name == "texture_normal")
                number = std::to_string(normalNr++); // transfer unsigned int to stream
            else if(name == "texture_height")
                number = std::to_string(heightNr++); // transfer unsigned int to stream
            else if(name == "texture_reflection")// We'll now also need to add the code to set and bind to reflection textures
                number = std::to_string(reflectionNr++);

            // now set the sampler to the correct texture unit
            glUniform1i(glGetUniformLocation(shader.ID, (name + number).c_str()), i);
            // and finally bind the texture
            glBindTexture(GL_TEXTURE_2D, textures[i].id);
        }
        
        // draw mesh
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);//draw完之后解绑防止误用

        // always good practice to set everything back to defaults once configured.
        glActiveTexture(GL_TEXTURE0);
    }

    bool isOnBoundary(HalfedgeHandle const & heh) {
        assert(heh.face_id != IndexNotValid);
        return heh.oppo_halfedge_id == IndexNotValid;// || heh.face_id == IndexNotValid;
    }

    bool isOnBoundary(VertexHandle const & vh) {
        for (unsigned i = 0; i < vh.outgoing_halfedge_ids.size(); i++){
            if (isOnBoundary(halfedgeList[vh.outgoing_halfedge_ids[i]]))
                return true;
        }
        for (unsigned i = 0; i < vh.incoming_halfedge_ids.size(); i++){
            if (isOnBoundary(halfedgeList[vh.incoming_halfedge_ids[i]]))
                return true;
        }
        return false;
    }

    bool is_collapse_OK(HalfedgeHandle const & v0v1) {
        // is the edge already deleted?
        if (edgeList[v0v1.edge_id].deleted)
            return false;

        // HalfedgeHandle const & v1v0 = halfedgeList[v0v1.oppo_halfedge_id];
        size_t v0_id = v0v1.from_vertex_id, v1_id = v0v1.to_vertex_id;

        // are vertices already deleted ?
        if (vertexList[v0_id].deleted || vertexList[v1_id].deleted)
            return false;

        size_t vl_id, vr_id;
        // the edges v1-vl and vl-v0 must not be both boundary edges
        if (!isOnBoundary(v0v1)) {//v0v1不在边界上则说明v1v0不为IndexNotValid
            HalfedgeHandle const & h1 = halfedgeList[v0v1.next_halfedge_id];
            HalfedgeHandle const & h2 = halfedgeList[h1.next_halfedge_id];
            vl_id = h1.to_vertex_id;
            // if (isOnBoundary(halfedgeList[h1.oppo_halfedge_id]) &&
            //     isOnBoundary(halfedgeList[h2.oppo_halfedge_id])){
            if (h1.oppo_halfedge_id == IndexNotValid && h2.oppo_halfedge_id == IndexNotValid){
                return false;
            }
        }

        // the edges v0-vr and vr-v1 must not be both boundary edges
        if (v0v1.oppo_halfedge_id != IndexNotValid) {
            HalfedgeHandle const & v1v0 = halfedgeList[v0v1.oppo_halfedge_id];
            HalfedgeHandle const & h1 = halfedgeList[v1v0.next_halfedge_id];
            HalfedgeHandle const & h2 = halfedgeList[h1.next_halfedge_id];
            vr_id = h1.to_vertex_id;
            // if (isOnBoundary(halfedgeList[h1.oppo_halfedge_id]) &&
            //     isOnBoundary(halfedgeList[h2.oppo_halfedge_id])){
            if (h1.oppo_halfedge_id == IndexNotValid && h2.oppo_halfedge_id == IndexNotValid){
                return false;
            }
        }

        // if vl and vr are equal or both invalid -> fail
        if (vl_id == vr_id) return false;

        // test intersection of the one-rings of v0 and v1
        std::unordered_set<size_t> v0_one_rings;
        VertexHandle const & v0 = vertexList[v0_id];
        VertexHandle const & v1 = vertexList[v1_id];
        for (unsigned i = 0; i < v0.outgoing_halfedge_ids.size(); i++)
            v0_one_rings.insert(halfedgeList[v0.outgoing_halfedge_ids[i]].to_vertex_id);
        for (unsigned i = 0; i < v0.incoming_halfedge_ids.size(); i++)
            v0_one_rings.insert(halfedgeList[v0.incoming_halfedge_ids[i]].from_vertex_id);

        std::unordered_set<size_t> v0v1_one_rings_intersect;
        for (unsigned i = 0; i < v1.outgoing_halfedge_ids.size(); i++){
            if (v0_one_rings.find(halfedgeList[v1.outgoing_halfedge_ids[i]].to_vertex_id) != v0_one_rings.end())
                v0v1_one_rings_intersect.insert(halfedgeList[v1.outgoing_halfedge_ids[i]].to_vertex_id);
        }
        for (unsigned i = 0; i < v1.incoming_halfedge_ids.size(); i++){
            if (v0_one_rings.find(halfedgeList[v1.incoming_halfedge_ids[i]].from_vertex_id) != v0_one_rings.end())
                v0v1_one_rings_intersect.insert(halfedgeList[v1.incoming_halfedge_ids[i]].from_vertex_id);
        }

        for (std::unordered_set<size_t>::iterator it = v0v1_one_rings_intersect.begin(); it != v0v1_one_rings_intersect.end(); it++){
            if (*it != vl_id && *it != vr_id)
                return false;
        }

        bool isOnBoundary_v1v0 = (v0v1.oppo_halfedge_id==IndexNotValid) ? true : isOnBoundary(halfedgeList[v0v1.oppo_halfedge_id]);
        // edge between two boundary vertices should be a boundary edge
        if ( isOnBoundary(v0)   &&  isOnBoundary(v1) &&
            !isOnBoundary(v0v1) && !isOnBoundary_v1v0)
            return false;

        // passed all tests
        return true;
    }

    void delete_halfedge(HalfedgeHandle & hh) {
        hh.deleted = true;
        // 更新和它有关系的顶点
        VertexHandle & v_to = vertexList[hh.to_vertex_id];
        assert(remove_from_vector(v_to.incoming_halfedge_ids, hh.halfedge_id));

        VertexHandle & v_from = vertexList[hh.from_vertex_id];
        assert(remove_from_vector(v_from.outgoing_halfedge_ids, hh.halfedge_id));

        // 更新和它有关系的所有半边！
        if (hh.next_halfedge_id != IndexNotValid){
            halfedgeList[hh.next_halfedge_id].prev_halfedge_id = IndexNotValid;
            hh.next_halfedge_id = IndexNotValid;
        }
        if (hh.prev_halfedge_id != IndexNotValid){
            halfedgeList[hh.prev_halfedge_id].next_halfedge_id = IndexNotValid;
            hh.prev_halfedge_id = IndexNotValid;
        }
        if (hh.oppo_halfedge_id != IndexNotValid){ 
            halfedgeList[hh.oppo_halfedge_id].oppo_halfedge_id = IndexNotValid;
            hh.oppo_halfedge_id = IndexNotValid;
        }
        // 更新和它有关系的边
        EdgeHandle & edge = edgeList[hh.edge_id];
        if (edge.halfedge_1_id == hh.halfedge_id)
            edge.halfedge_1_id = IndexNotValid;
        else if (edge.halfedge_2_id == hh.halfedge_id)
            edge.halfedge_2_id = IndexNotValid;
        else {
            std::cout << "Error: could not find corresponding edge when deleting halfedge " << hh.halfedge_id << "!" << std::endl << std::endl;
            exit(1);
        }
        if (edge.halfedge_1_id == IndexNotValid && edge.halfedge_2_id == IndexNotValid){//这条边已经没有了存在的意义
            edge.deleted = true;
            hh.edge_id = IndexNotValid;
        }
    }

    void make_halfedges_pair(HalfedgeHandle & host, HalfedgeHandle & guest) {
        assert(host.oppo_halfedge_id == IndexNotValid);
        host.oppo_halfedge_id = guest.halfedge_id;

        assert(guest.oppo_halfedge_id == IndexNotValid);
        guest.oppo_halfedge_id = host.halfedge_id;

        guest.edge_id = host.edge_id;
        // 更新边的信息
        if (edgeList[host.edge_id].halfedge_1_id == host.halfedge_id)
            edgeList[host.edge_id].halfedge_2_id = guest.halfedge_id;
        else if (edgeList[host.edge_id].halfedge_2_id == host.halfedge_id)
            edgeList[host.edge_id].halfedge_1_id = guest.halfedge_id;
        assert(edgeList[host.edge_id].halfedge_1_id != IndexNotValid && edgeList[host.edge_id].halfedge_2_id != IndexNotValid);
        // 修改guest的前后半边信息不用修改！
    }

    void delete_around(HalfedgeHandle & hh) {
        size_t start_he_id = hh.halfedge_id, he_id = start_he_id;
        do {
            HalfedgeHandle & he_del = halfedgeList[he_id];
            he_id = he_del.next_halfedge_id;
            delete_halfedge(he_del);
        } while (he_id != IndexNotValid);
    }

    void zyCollapse(HalfedgeHandle & _hh) {
        // remove faces
        size_t oppo_he_id = _hh.oppo_halfedge_id;
        delete_around(_hh);
        if (oppo_he_id != IndexNotValid) {
            delete_around(halfedgeList[oppo_he_id]);
        }
        // move halfedges
        VertexHandle & v_to   = vertexList[_hh.to_vertex_id];
        VertexHandle & v_from = vertexList[_hh.from_vertex_id]; v_from.deleted = true;
        // 将原指向v_from的改为指向v_to
        for (unsigned i = 0; i < v_from.incoming_halfedge_ids.size(); i++){
            HalfedgeHandle & he_p2_vfrom = halfedgeList[v_from.incoming_halfedge_ids[i]];// halfedge point to v_from
            if (!he_p2_vfrom.deleted){//这条边仍有效
                he_p2_vfrom.to_vertex_id = v_to.vertex_id;
                v_to.incoming_halfedge_ids.push_back(he_p2_vfrom.halfedge_id);
                if (he_p2_vfrom.oppo_halfedge_id == IndexNotValid) {
                    // 如果原指向v_from的这条半边已经没有了对半边，则从v_to的出射半边中找，而边的归属归新的边
                    for (unsigned i = 0; i < v_to.outgoing_halfedge_ids.size(); i++){// 找到和它重新匹配的半边
                        HalfedgeHandle & he_pf_vto = halfedgeList[v_to.outgoing_halfedge_ids[i]];//halfedge point from v_to
                        if (he_pf_vto.to_vertex_id == he_p2_vfrom.from_vertex_id) 
                            make_halfedges_pair(he_pf_vto, he_p2_vfrom);
                    }
                } else {// 如果原指向v_from的这条半边仍然有对半边，则它和它的对半边不动，全部移交给v_to

                }
            }
        }
        // 将原从v_from指出的改为从v_tp指出
        for (unsigned i = 0; i < v_from.outgoing_halfedge_ids.size(); i++){
            HalfedgeHandle & he_pf_vfrom = halfedgeList[v_from.outgoing_halfedge_ids[i]];// halfedge point from v_from
            if (!he_pf_vfrom.deleted){//这条边仍有效
                he_pf_vfrom.from_vertex_id = v_to.vertex_id;
                v_to.outgoing_halfedge_ids.push_back(he_pf_vfrom.halfedge_id);
                if (he_pf_vfrom.oppo_halfedge_id == IndexNotValid){
                    // 如果原从v_from指出的这条半边已经没有了对半边，则从v_to的入射半边中找，而边的归属归新的边
                    for (unsigned i = 0; i < v_to.incoming_halfedge_ids.size(); i++){// 找到和它重新匹配的半边
                        HalfedgeHandle & he_p2_vto = halfedgeList[v_to.incoming_halfedge_ids[i]];// halfedge point to v_to
                        if (he_p2_vto.from_vertex_id == he_pf_vfrom.to_vertex_id)
                            make_halfedges_pair(he_p2_vto, he_pf_vfrom);
                    }
                } else {

                }
            }
        }
        

    }

private:
    bool remove_from_vector(std::vector<size_t> & vec, size_t const & target) {
        for (std::vector<size_t>::iterator it = vec.begin(); it != vec.end(); it++){
            if (*it == target){
                vec.erase(it);
                return true;
            }
        }
        return false;
    }

    void loadData(std::string filename) {
        std::ifstream fin;
        fin.open(filename.c_str(), std::ios::in);
        std::string sline;
        char s0;
        vec3 point;
        std::vector<size_t> vids(3, IndexNotValid);
        
        while (std::getline(fin, sline)) {//从指定文件逐行读取
            if (sline[0] == 'v') {
                std::istringstream ins(sline);
                ins >> s0 >> point.x >> point.y >> point.z;
                // sscanf(sline.c_str(), "%c %f %f %f", &s0, &point.x, &point.y, &point.z);
                // ins >> s0 >> point.x >> point.y >> point.z;
                add_vertex(point);
            } else if (sline[0] == 'f') {
                std::istringstream ins(sline);
                ins >> s0 >> vids[0] >> vids[1] >> vids[2];
                // scanf(sline.c_str(), "%c %d %d %d", &s0, &vids[0] , &vids[1] , &vids[2]);
                vids[0] -= 1; vids[1] -= 1; vids[2] -= 1;
                add_face(vids);
            }
        }
        fin.close();
    }
};

#endif