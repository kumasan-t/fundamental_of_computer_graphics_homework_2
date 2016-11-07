#include "tesselation.h"

// make normals for each face - duplicates all vertex data
void facet_normals(Mesh* mesh) {
    // allocates new arrays
    auto pos = vector<vec3f>();
    auto norm = vector<vec3f>();
    auto texcoord = vector<vec2f>();
    auto triangle = vector<vec3i>();
    auto quad = vector<vec4i>();
    // froeach triangle
    for(auto f : mesh->triangle) {
        // grab current pos size
        auto nv = (int)pos.size();
        // compute face face normal
        auto fn = normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x]));
        // add triangle
        triangle.push_back({nv,nv+1,nv+2});
        // add vertex data
        for(auto i : range(3)) {
            pos.push_back(mesh->pos[f[i]]);
            norm.push_back(fn);
            if(not mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
        }
    }
    // froeach quad
    for(auto f : mesh->quad) {
        // grab current pos size
        auto nv = (int)pos.size();
        // compute face normal
        auto fn = normalize(normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x])) +
                            normalize(cross(mesh->pos[f.z]-mesh->pos[f.x], mesh->pos[f.w]-mesh->pos[f.x])));
        // add quad
        quad.push_back({nv,nv+1,nv+2,nv+3});
        // add vertex data
        for(auto i : range(4)) {
            pos.push_back(mesh->pos[f[i]]);
            norm.push_back(fn);
            if(not mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
        }
    }
    // set back mesh data
    mesh->pos = pos;
    mesh->norm = norm;
    mesh->texcoord = texcoord;
    mesh->triangle = triangle;
    mesh->quad = quad;
}

// smooth out normal - does not duplicate data
void smooth_normals(Mesh* mesh) {
    // PLACEHOLDER CODE - REMOVE AFTER FUNCTION IS IMPLEMENTED
    // YOUR CODE GOES HERE ---------------------
    // set normals array to the same length as pos and init all elements to zero
    mesh->norm.resize(mesh->pos.size());
    // foreach triangle
	vector<vec3f> normals = mesh->norm;
	for (vec3i triangle : mesh->triangle) {
		// compute face normal
		vec3f first_side_vector = mesh->pos[triangle.x] - mesh->pos[triangle.y];
		vec3f second_side_vector = mesh->pos[triangle.x] - mesh->pos[triangle.z];
		vec3f triangle_face_normal = cross(first_side_vector, second_side_vector);
		// accumulate face normal to the vertex normals of each face index
		normals[triangle.x] += triangle_face_normal;
		normals[triangle.y] += triangle_face_normal;
		normals[triangle.z] += triangle_face_normal;
	}
	// foreach quad
	for (vec4i quad : mesh->quad) {
		// compute face normal
		vec3f quad_face_normal = normalize(
			normalize(cross(mesh->pos[quad.y] - mesh->pos[quad.x], mesh->pos[quad.z] - mesh->pos[quad.x])) +
			normalize(cross(mesh->pos[quad.z] - mesh->pos[quad.x], mesh->pos[quad.w] - mesh->pos[quad.x])));
		// accumulate face normal to the vertex normals of each face index
		normals[quad.x] += quad_face_normal;
		normals[quad.y] += quad_face_normal;
		normals[quad.w] += quad_face_normal;
		normals[quad.z] += quad_face_normal;
	}
    // normalize all vertex normals
	for (int i = 0; i < normals.size(); i++) {
		normals[i] = normalize(normals[i]);
	}
	mesh->norm = normals;
}

// smooth out tangents
void smooth_tangents(Mesh* polyline) {
    // set tangent array
    polyline->norm = vector<vec3f>(polyline->pos.size(),zero3f);
    // foreach line
    for(auto l : polyline->line) {
        // compute line tangent
        auto lt = normalize(polyline->pos[l.y]-polyline->pos[l.x]);
        // accumulate segment tangent to vertex tangent on each vertex
        for (auto i : range(2)) polyline->norm[l[i]] += lt;
    }
    // normalize all vertex tangents
    for (auto& t : polyline->norm) t = normalize(t);
}

// apply Catmull-Clark mesh subdivision
// does not subdivide texcoord
void subdivide_catmullclark(Mesh* subdiv) {
    // YOUR CODE GOES HERE ---------------------
    // skip is needed
    // allocate a working Mesh copied from the subdiv
	Mesh* working_mesh = subdiv;
    // foreach level
	for (int i = 0; i < working_mesh->subdivision_catmullclark_level; i++) {
		// make empty pos and quad arrays
		vector<vec3f> working_pos;
		vector<vec4i> working_quad;
		// create edge_map from current mesh
		EdgeMap edge_map = EdgeMap(working_mesh->triangle,working_mesh->quad);
		// linear subdivision - create vertices
		vector<vec3f> working_vertices;
		// copy all vertices from the current mesh
		for (auto triangle : working_mesh->triangle) {
			working_vertices.push_back(working_mesh->pos[triangle.x]);
			working_vertices.push_back(working_mesh->pos[triangle.y]);
			working_vertices.push_back(working_mesh->pos[triangle.z]);
		}
		for (auto quad : working_mesh->quad) {
			working_vertices.push_back(working_mesh->pos[quad.x]);
			working_vertices.push_back(working_mesh->pos[quad.y]);
			working_vertices.push_back(working_mesh->pos[quad.w]);
			working_vertices.push_back(working_mesh->pos[quad.z]);
		}
		// add vertices in the middle of each edge (use EdgeMap)
		for (vec2i edge : edge_map._edge_list) {
			vec3f temp_mid_edge_vertex = (working_mesh->pos[edge.x] + working_mesh->pos[edge.y]) / 2.0f;
			working_pos.push_back(temp_mid_edge_vertex);
		}
		// add vertices in the middle of each triangle
		for (auto triangle : working_mesh->triangle) {
			vec3f temp_mid_triangle_vertex = (
				working_mesh->pos[triangle.x] +
				working_mesh->pos[triangle.y] +
				working_mesh->pos[triangle.z])
				/ 3.0f;
			working_pos.push_back(temp_mid_triangle_vertex);
		}
		// add vertices in the middle of each quad
		for (auto quad : working_mesh->quad) {
			vec3f temp_mid_quad_vertex = (
				working_mesh->pos[quad.x] +
				working_mesh->pos[quad.y] +
				working_mesh->pos[quad.w] +
				working_mesh->pos[quad.z] )
				/ 4.0f;
			working_pos.push_back(temp_mid_quad_vertex);
		}
		// subdivision pass --------------------------------
		// compute an offset for the edge vertices
		int edge_vertex_index = working_mesh->triangle.size() + working_mesh->quad.size() - 1;
		// compute an offset for the triangle vertices
		int triangle_vertex_index = edge_vertex_index + edge_map._edge_list.size();
		// compute an offset for the quad vertices
		int quad_vertex_index = triangle_vertex_index + working_mesh->triangle.size();
		// foreach triangle
		for (auto triangle : working_mesh->triangle) {
			// add three quads to the new quad array	
			vec4i temp_quad = vec4i(
				triangle.x,
				edge_vertex_index,
				triangle_vertex_index,
				edge_vertex_index + 2
				);
			vec4i temp_quad = vec4i(
				triangle.y,
				edge_vertex_index,
				triangle_vertex_index,
				edge_vertex_index + 1
				);
			vec4i temp_quad = vec4i(
				triangle.z,
				edge_vertex_index + 2,
				triangle_vertex_index,
				edge_vertex_index + 1
				);

			edge_vertex_index += 3;
			triangle_vertex_index += 1;
			working_quad.push_back(temp_quad);
		}
		// foreach quad
		// add four quads to the new quad array
		// averaging pass ----------------------------------
		// create arrays to compute pos averages (avg_pos, avg_count)
		// arrays have the same length as the new pos array, and are init to zero
		// for each new quad
		// compute quad center using the new pos array
		// foreach vertex index in the quad
		// normalize avg_pos with its count avg_count
		// correction pass ----------------------------------
		// foreach pos, compute correction p = p + (avg_p - p) * (4/avg_count)
		// set new arrays pos, quad back into the working mesh; clear triangle array
	}
    // clear subdivision
    // according to smooth, either smooth_normals or facet_normals
    // copy back
    // clear
}

// subdivide bezier spline into line segments (assume bezier has only bezier segments and no lines)
void subdivide_bezier(Mesh* bezier) {
    // YOUR CODE GOES HERE ---------------------
    // skip is needed
    // allocate a working polyline from bezier
    // foreach level
        // make new arrays of positions and bezier segments
        // copy all the vertices into the new array (this waste space but it is easier for now)
        // foreach bezier segment
            // apply subdivision algorithm
            // prepare indices for two new segments
            // add mid point
            // add points for first segment and fix segment indices
            // add points for second segment and fix segment indices
            // add indices for both segments into new segments array
        // set new arrays pos, segments into the working lineset
    // copy bezier segments into line segments
    // clear bezier array from lines
    // run smoothing to get proper tangents
    // copy back
    // clear
}

Mesh* make_surface_mesh(frame3f frame, float radius, bool isquad, Material* mat, float offset) {
    auto mesh = new Mesh{};
    mesh->frame = frame;
    mesh->mat = mat;
    if(isquad) {
        mesh->pos = { {-radius,-radius,-offset}, {radius,-radius,-offset},
            {radius,radius,-offset}, {-radius,radius,-offset} };
        mesh->norm = {z3f,z3f,z3f,z3f};
        mesh->quad = { {0,1,2,3} };
    } else {
        map<pair<int,int>,int> vid;
        for(auto j : range(64+1)) {
            for(auto i : range(128+1)) {
                auto u = 2 * pif * i / 64.0f, v = pif * j / 32.0f;
                auto d = vec3f{cos(u)*sin(v),sin(u)*sin(v),cos(v)};
                vid[{i,j}] = mesh->pos.size();
                mesh->pos.push_back(d*radius*(1-offset));
                mesh->norm.push_back(d);
            }
        }
        for(auto j : range(64)) {
            for(auto i : range(128)) {
                mesh->quad.push_back({vid[{i,j}],vid[{i+1,j}],vid[{i+1,j+1}],vid[{i,j+1}]});
            }
        }
    }
    return mesh;
}

void subdivide_surface(Surface* surface) {
    surface->_display_mesh = make_surface_mesh(
        surface->frame, surface->radius, surface->isquad, surface->mat);
}

void subdivide(Scene* scene) {
    for(auto mesh : scene->meshes) {
        if(mesh->subdivision_catmullclark_level) subdivide_catmullclark(mesh);
        if(mesh->subdivision_bezier_level) subdivide_bezier(mesh);
    }
    for(auto surface : scene->surfaces) {
        subdivide_surface(surface);
    }
}
