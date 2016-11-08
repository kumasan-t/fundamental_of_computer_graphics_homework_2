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
		vector<vec3f> pos;
		vector<vec4i> new_quad;
		// create edge_map from current mesh
		EdgeMap edge_map = EdgeMap(working_mesh->triangle,working_mesh->quad);
		// linear subdivision - create vertices
		// copy all vertices from the current mesh
		pos = working_mesh->pos;
		// add vertices in the middle of each edge (use EdgeMap)
		for (vec2i edge : edge_map.edges()) {
			vec3f temp_mid_edge_vertex = (working_mesh->pos[edge.x] + working_mesh->pos[edge.y]) / 2.0f;
			pos.push_back(temp_mid_edge_vertex);
		}
		// add vertices in the middle of each triangle
		for (vec3i triangle : working_mesh->triangle) {
			vec3f temp_mid_triangle_vertex = (
				working_mesh->pos[triangle.x] +
				working_mesh->pos[triangle.y] +
				working_mesh->pos[triangle.z])
				/ 3.0f;
			pos.push_back(temp_mid_triangle_vertex);
		}
		// add vertices in the middle of each quad
		for (vec4i quad : working_mesh->quad) {
			vec3f temp_mid_quad_vertex = (
				working_mesh->pos[quad.x] +
				working_mesh->pos[quad.y] +
				working_mesh->pos[quad.z] +
				working_mesh->pos[quad.w] )
				/ 4.0f;
			pos.push_back(temp_mid_quad_vertex);
		}
		// subdivision pass --------------------------------
		// compute an offset for the edge vertices
		int edge_vertex_index = working_mesh->pos.size();
		// compute an offset for the triangle vertices
		int triangle_vertex_index = edge_vertex_index + edge_map.edges().size();
		// compute an offset for the quad vertices
		int quad_vertex_index = triangle_vertex_index + working_mesh->triangle.size();
		// foreach triangle
		for (vec3i triangle : working_mesh->triangle) {
			// add three quads to the new quad array	
			new_quad.push_back(vec4i(
				triangle.x,
				edge_vertex_index + edge_map.edge_index(vec2i(triangle.x,triangle.y)),
				triangle_vertex_index,
				edge_vertex_index + edge_map.edge_index(vec2i(triangle.x, triangle.z))
				));
			new_quad.push_back(vec4i(
				triangle.y,
				edge_vertex_index + edge_map.edge_index(vec2i(triangle.y,triangle.z)),
				triangle_vertex_index,
				edge_vertex_index + edge_map.edge_index(vec2i(triangle.y, triangle.x))
				));
			new_quad.push_back(vec4i(
				triangle.z,
				edge_vertex_index + edge_map.edge_index(vec2i(triangle.z, triangle.x)),
				triangle_vertex_index,
				edge_vertex_index + edge_map.edge_index(vec2i(triangle.z, triangle.y))
				));
			triangle_vertex_index++;
		}
		// foreach quad
		for (vec4i quad : working_mesh->quad) {
			// add four quads to the new quad array
			new_quad.push_back(vec4i(
				quad.w,
				edge_vertex_index + edge_map.edge_index(vec2i(quad.w, quad.x)),
				quad_vertex_index,
				edge_vertex_index + edge_map.edge_index(vec2i(quad.w, quad.z))
				));
			new_quad.push_back(vec4i(
				quad.x,
				edge_vertex_index + edge_map.edge_index(vec2i(quad.x, quad.y)),
				quad_vertex_index,
				edge_vertex_index + edge_map.edge_index(vec2i(quad.x, quad.w))
				));
			new_quad.push_back(vec4i(
				quad.y,
				edge_vertex_index + edge_map.edge_index(vec2i(quad.y, quad.z)),
				quad_vertex_index,
				edge_vertex_index + edge_map.edge_index(vec2i(quad.y, quad.x))
				));
			new_quad.push_back(vec4i(
				quad.z,
				edge_vertex_index + edge_map.edge_index(vec2i(quad.z, quad.w)),
				quad_vertex_index,
				edge_vertex_index + edge_map.edge_index(vec2i(quad.z, quad.y))
				));
			quad_vertex_index++;
		}
		// averaging pass ----------------------------------
		// create arrays to compute pos averages (avg_pos, avg_count)
		vector<vec3f> avg_pos;
		vector<int> avg_count;
		// arrays have the same length as the new pos array, and are init to zero
		avg_pos.resize(pos.size());
		avg_count.resize(pos.size());
		for (int i = 0; i < avg_pos.size(); i++) { avg_pos[i] = zero3f; }
		for (int i = 0; i < avg_count.size(); i++) { avg_count[i] = 0;  }
		// for each new quad
		for (vec4i quad : new_quad) {
			// compute quad center using the new pos array
			vec3f temp_mid_new_quad_vertex = (
				pos[quad.x] +
				pos[quad.y] +
				pos[quad.z] +
				pos[quad.w])
				/ 4.0f;
			// foreach vertex index in the quad
			for (int i = 0; i < 4; i++) {
				avg_pos[quad[i]] += temp_mid_new_quad_vertex;
				avg_count[quad[i]] += 1;
			}
		}
		for (int i : range(avg_count.size())) {
			// normalize avg_pos with its count avg_count
			avg_pos[i] /= (float)avg_count[i];
		}
		// correction pass ----------------------------------
		// foreach pos, compute correction p = p + (avg_p - p) * (4/avg_count)
		for (int i : range(pos.size())) {
			pos[i] += (avg_pos[i] - pos[i]) * (4.0f / avg_count[i]);
		}
		// set new arrays pos, quad back into the working mesh; clear triangle array
		working_mesh->pos = pos;
		working_mesh->quad = new_quad;
		working_mesh->triangle.clear();
		message("ARRIVATO FIN QUI!!!1\n");
	}
    // clear subdivision
    // according to smooth, either smooth_normals or facet_normals
	message("ARRIVATO FIN QUI!!!\n");
	(working_mesh->subdivision_catmullclark_smooth) ? smooth_normals(working_mesh) : facet_normals(working_mesh);
    // copy back
	subdiv = working_mesh;
    // clear
}

// subdivide bezier spline into line segments (assume bezier has only bezier segments and no lines)
void subdivide_bezier(Mesh* bezier) {
    // YOUR CODE GOES HERE ---------------------
    // skip is needed
    // allocate a working polyline from bezier
	Mesh* working_polyline = bezier;
    // foreach level
	for (int i = 0; i < bezier->subdivision_bezier_level; i++) {
        // make new arrays of positions and bezier segments
		vector<vec3f> pos;
		vector<vec4i> new_spline;
        // copy all the vertices into the new array (this waste space but it is easier for now)
		pos = working_polyline->pos;
        // foreach bezier segment
		for (vec4i b_spline : working_polyline->spline) {
			// apply subdivision algorithm
			auto r_0 = (pos[b_spline[0]] + pos[b_spline[1]]) / 2;
			auto r_1 = (pos[b_spline[1]] + pos[b_spline[2]]) / 2;
			auto r_2 = (pos[b_spline[2]] + pos[b_spline[3]]) / 2;

			auto s_0 = (r_0 + r_1) / 2;
			auto s_1 = (r_1 + r_2) / 2;

			auto m_0 = (s_0 + s_1) / 2;
 			// prepare indices for two new segments
			auto r_0_index = pos.size();
			pos.push_back(r_0);

			auto r_2_index = pos.size();
			pos.push_back(r_2);

			// add mid point
			auto s_0_index = pos.size();
			pos.push_back(s_0);

			auto s_1_index = pos.size();
			pos.push_back(s_1);

			auto m_0_index = pos.size();
			pos.push_back(m_0);
			// add points for first segment and fix segment indices
			vec4i subdiv_spline_1 = vec4i(b_spline[0], r_0_index, s_0_index, m_0_index);
			// add points for second segment and fix segment indices
			vec4i subdiv_spline_2 = vec4i(m_0_index, s_1_index, r_2_index, b_spline[3]);
			// add indices for both segments into new segments array
			new_spline.push_back(subdiv_spline_1);
			new_spline.push_back(subdiv_spline_2);
		}
        // set new arrays pos, segments into the working lineset
		working_polyline->pos = pos;
		working_polyline->spline = new_spline;
	}
    // copy bezier segments into line segments
	for (vec4i spline : working_polyline->spline) {
		working_polyline->line.push_back(vec2i(spline[0], spline[1]));
		working_polyline->line.push_back(vec2i(spline[1], spline[2]));
		working_polyline->line.push_back(vec2i(spline[2], spline[3]));
	}
    // clear bezier array from lines
	bezier->line.clear();
    // run smoothing to get proper tangents
	smooth_tangents(working_polyline);
    // copy back
	bezier = working_polyline;
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
