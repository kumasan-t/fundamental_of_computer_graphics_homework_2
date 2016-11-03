#version 120

varying vec3 pos;                   // [from vertex shader] position in world space
varying vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
varying vec2 texcoord;              // [from vertex shader] texture coordinate

uniform vec3 camera_pos;            // camera position (center of the camera frame)

uniform vec3 ambient;               // scene ambient

uniform int lights_num;             // number of lights
uniform vec3 light_pos[16];         // light positions
uniform vec3 light_intensity[16];   // light intensities

uniform vec3 material_kd;           // material kd
uniform vec3 material_ks;           // material ks
uniform float material_n;           // material n
uniform bool material_is_lines;     // whether the material is lines or meshes

uniform bool material_kd_txt_on;    // material kd texture enabled
uniform sampler2D material_kd_txt;  // material kd texture
uniform bool material_ks_txt_on;    // material ks texture enabled
uniform sampler2D material_ks_txt;  // material ks texture
uniform bool material_norm_txt_on;    // material norm texture enabled
uniform sampler2D material_norm_txt;  // material norm texture

float distSqr(vec3 a, vec3 b) { return pow(length(a - b),2); }
// main
void main() {
	vec3 ks = material_ks;
	vec3 kd = material_kd;
    // re-normalize normals
    vec3 n = normalize(norm);
    vec3 accumulated_color = vec3(1,0,0);   // initialize to red to see it well
    // YOUR CODE GOES HERE ---------------------
    // lookup normal map if needed
    // compute material values by looking up textures is necessary
	if (material_kd_txt_on) { kd = material_kd * vec3(texture2D(material_kd_txt,texcoord)); }
	if (material_ks_txt_on) { ks = material_ks * vec3(texture2D(material_ks_txt,texcoord)); }
	if (material_norm_txt_on) { n = normalize(2 * vec3(texture2D(material_norm_txt,texcoord)) - vec3(1)); }
    // accumulate ambient
	accumulated_color = ambient * kd;
    // foreach light
	for (int i = 0; i < lights_num; i++) {
        // compute point light color at pos
		vec3 light_color = light_intensity[i] / distSqr(light_pos[i],pos);
        // compute light direction at pos
		vec3 light_direction = normalize(light_pos[i] - pos);
        // compute view direction using camera_pos and pos
		vec3 viewer_direction = normalize(camera_pos - pos);
        // compute h
		vec3 bisector_h = normalize(light_direction + viewer_direction);
		vec3 diffuse = kd;
		vec3 specular = ks * pow(max(0.0f, dot(bisector_h,n)), material_n);
        // accumulate blinn-phong model
		accumulated_color += light_color * (diffuse + specular) * max(0,dot(n,light_direction));
	}
    // output final color by setting gl_FragColor
    gl_FragColor = vec4(accumulated_color,1);
}
