#include "ViewerApplication.hpp"

#include <iostream>
#include <numeric>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/io.hpp>

#include "utils/CameraController.hpp"
#include "utils/cameras.hpp"
#include "utils/gltf.hpp"
#include "utils/images.hpp"

#include <stb_image_write.h>
#include <tiny_gltf.h>

const GLuint VERTEX_ATTRIB_POSITION_IDX = 0;
const GLuint VERTEX_ATTRIB_NORMAL_IDX = 1;
const GLuint VERTEX_ATTRIB_TEXCOORD0_IDX = 2;
const GLuint VERTEX_ATTRIB_TANGENT_IDX = 3;

// Used only to clean up the buffers.
std::vector<GLuint> tangentBuffers;

void keyCallback(
    GLFWwindow *window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE) {
    glfwSetWindowShouldClose(window, 1);
  }
}

int ViewerApplication::run()
{
  // Loader shaders
  const auto glslProgram = compileProgram({m_ShadersRootPath / m_vertexShader,
      m_ShadersRootPath / m_fragmentShader});

  const auto modelViewProjMatrixLocation =
      glGetUniformLocation(glslProgram.glId(), "uModelViewProjMatrix");
  const auto modelViewMatrixLocation =
      glGetUniformLocation(glslProgram.glId(), "uModelViewMatrix");
  const auto normalMatrixLocation =
      glGetUniformLocation(glslProgram.glId(), "uNormalMatrix");

  const auto lightDirectionLocation =
      glGetUniformLocation(glslProgram.glId(), "uLightDirection");
  const auto lightIntensityLocation =
      glGetUniformLocation(glslProgram.glId(), "uLightIntensity");

  const auto baseColorTextureLocation =
      glGetUniformLocation(glslProgram.glId(), "uBaseColorTexture");
  const auto baseColorFactorLocation =
      glGetUniformLocation(glslProgram.glId(), "uBaseColorFactor");

  const auto metallicRoughnessTextureLocation =
      glGetUniformLocation(glslProgram.glId(), "uMetallicRoughnessTexture");
  const auto metallicFactorLocation =
      glGetUniformLocation(glslProgram.glId(), "uMetallicFactor");
  const auto roughnessFactorLocation =
      glGetUniformLocation(glslProgram.glId(), "uRoughnessFactor");

  const auto emissiveTextureLocation =
      glGetUniformLocation(glslProgram.glId(), "uEmissiveTexture");
  const auto emissiveFactorLocation =
      glGetUniformLocation(glslProgram.glId(), "uEmissiveFactor");

  const auto occlusionTextureLocation =
      glGetUniformLocation(glslProgram.glId(), "uOcclusionTexture");
  const auto occlusionFactorLocation =
      glGetUniformLocation(glslProgram.glId(), "uOcclusionFactor");
  const auto occlusionEnabledLocation =
      glGetUniformLocation(glslProgram.glId(), "uOcclusionEnabled");

  const auto normalMapTextureLocation =
      glGetUniformLocation(glslProgram.glId(), "uNormalMapTexture");

  glm::vec3 lightDirection(glm::sin(0.f) * glm::cos(0.f), glm::cos(0.f),
      glm::sin(0.f) * glm::sin(0.f));
  glm::vec3 lightIntensity(1.f, 1.f, 1.f);
  bool ambiantOcclusion = true;

  // Load the scene.
  tinygltf::Model model;
  // Loading the glTF file
  if (!loadGltfFile(model))
    throw std::exception("Cannot load the glTF model");

  // Creation of Texture Objects
  auto textureObjects = createTextureObjects(model);

  // Default texture.
  float white[] = {1, 1, 1, 1};
  GLuint whiteTexture = 0;
  glGenTextures(1, &whiteTexture);

  glBindTexture(GL_TEXTURE_2D, whiteTexture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, 1, 0, GL_RGBA, GL_FLOAT, white);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);

  glBindTexture(GL_TEXTURE_2D, 0);

  // Creation of Buffer Objects
  auto bufferObjects = createBufferObjects(model);

  // Creation of Vertex Array Objects
  std::vector<VaoRange> meshIndexToVaoRange;
  auto vertexArrayObjects =
      createVertexArrayObjects(model, bufferObjects, meshIndexToVaoRange);

  // Compute Scene BBox.
  glm::vec3 bboxMin, bboxMax;
  computeSceneBounds(model, bboxMin, bboxMax);
  glm::vec3 diag = (bboxMax - bboxMin);

  // Build projection matrix
  auto maxDistance = glm::length(diag);
  maxDistance = maxDistance > 0.f ? maxDistance : 100.f;
  const auto projMatrix =
      glm::perspective(70.f, float(m_nWindowWidth) / m_nWindowHeight,
          0.001f * maxDistance, 1.5f * maxDistance);

  std::unique_ptr<CameraController> cameraController =
      std::make_unique<TrackballCameraController>(
          m_GLFWHandle.window(), 0.6f * maxDistance);

  if (m_hasUserCamera) {
    cameraController->setCamera(m_userCamera);
  } else {
    glm::vec3 center = (bboxMin + bboxMax) / 2.f;
    glm::vec3 up = glm::vec3(0, 1, 0);
    glm::vec3 eye =
        (diag.z == 0.f) ? center + 2.f * glm::cross(diag, up) : center + diag;

    cameraController->setCamera(Camera{eye, center, up});
  }

  // Setup OpenGL state for rendering
  glEnable(GL_DEPTH_TEST);
  glslProgram.use();

  const auto bindMaterial = [&](const auto materialIndex) {
    if (materialIndex >= 0) {
      const auto &material = model.materials[materialIndex];
      const auto &pbrMetallicRoughness = material.pbrMetallicRoughness;

      // Base Color.
      const auto baseColorIndex = pbrMetallicRoughness.baseColorTexture.index;
      auto textureObject = whiteTexture;

      if (baseColorIndex >= 0) {
        textureObject = textureObjects[baseColorIndex];
      }

      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, textureObject);
      glUniform1i(baseColorTextureLocation, 0);

      glUniform4f(baseColorFactorLocation,
          (float)pbrMetallicRoughness.baseColorFactor[0],
          (float)pbrMetallicRoughness.baseColorFactor[1],
          (float)pbrMetallicRoughness.baseColorFactor[2],
          (float)pbrMetallicRoughness.baseColorFactor[3]);

      // Metallic/Roughness.
      textureObject = 0u;
      const auto baseMetallicRoughnessIndex =
          pbrMetallicRoughness.metallicRoughnessTexture.index;
      if (baseMetallicRoughnessIndex >= 0) {
        textureObject = textureObjects[baseMetallicRoughnessIndex];
      }

      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, textureObject);
      glUniform1i(metallicRoughnessTextureLocation, 1);

      glUniform1f(
          metallicFactorLocation, (float)pbrMetallicRoughness.metallicFactor);
      glUniform1f(
          roughnessFactorLocation, (float)pbrMetallicRoughness.roughnessFactor);

      // Emissive.
      textureObject = 0u;
      const auto emissiveTextureIndex = material.emissiveTexture.index;
      if (emissiveTextureIndex >= 0) {
        textureObject = textureObjects[emissiveTextureIndex];
      }

      glActiveTexture(GL_TEXTURE2);
      glBindTexture(GL_TEXTURE_2D, textureObject);
      glUniform1i(emissiveTextureLocation, 2);

      glUniform3f(emissiveFactorLocation, (float)material.emissiveFactor[0],
          (float)material.emissiveFactor[1], (float)material.emissiveFactor[2]);

      // Occlusion.
      textureObject = 0u;
      const auto occlusionTextureIndex = material.occlusionTexture.index;
      if (occlusionTextureIndex >= 0) {
        textureObject = textureObjects[occlusionTextureIndex];
      }

      glActiveTexture(GL_TEXTURE3);
      glBindTexture(GL_TEXTURE_2D, textureObject);
      glUniform1i(occlusionTextureLocation, 3);

      glUniform1f(
          occlusionFactorLocation, (float)material.occlusionTexture.strength);
      glUniform1i(occlusionEnabledLocation, (int)ambiantOcclusion);

      // Normal map.
      textureObject = 0u;
      const auto normalMapTextureIndex = material.normalTexture.index;
      if (normalMapTextureIndex >= 0) {
        textureObject = textureObjects[normalMapTextureIndex];
      }

      glActiveTexture(GL_TEXTURE4);
      glBindTexture(GL_TEXTURE_2D, textureObject);
      glUniform1i(normalMapTextureLocation, 4);

    } else {
      // Base color.
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, whiteTexture);
      glUniform1i(baseColorTextureLocation, 0);

      glUniform4f(baseColorFactorLocation, 1.f, 1.f, 1.f, 1.f);

      // Metallic/Roughness.
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, 0);
      glUniform1i(metallicRoughnessTextureLocation, 1);

      glUniform1f(metallicFactorLocation, 1.f);
      glUniform1f(roughnessFactorLocation, 1.f);

      // Emissive.
      glActiveTexture(GL_TEXTURE2);
      glBindTexture(GL_TEXTURE_2D, 0);
      glUniform1i(emissiveTextureLocation, 2);

      glUniform3f(emissiveFactorLocation, 0.f, 0.f, 0.f);

      // Occlusion.
      glActiveTexture(GL_TEXTURE3);
      glBindTexture(GL_TEXTURE_2D, 0);
      glUniform1i(occlusionTextureLocation, 3);

      glUniform1f(occlusionFactorLocation, 0.f);
      glUniform1i(occlusionEnabledLocation, 0);
    }
  };

  // Lambda function to draw the scene
  const auto drawScene = [&](const Camera &camera) {
    glViewport(0, 0, m_nWindowWidth, m_nWindowHeight);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    const auto viewMatrix = camera.getViewMatrix();

    if (lightDirectionLocation >= 0 && lightIntensityLocation >= 0) {
      const auto lightDirectionToSend =
          glm::normalize(glm::vec3(viewMatrix * glm::vec4(lightDirection, 0)));

      glUniform3f(lightIntensityLocation, lightIntensity.x, lightIntensity.y,
          lightIntensity.z);
      glUniform3f(lightDirectionLocation, lightDirectionToSend.x,
          lightDirectionToSend.y, lightDirectionToSend.z);
    }

    // The recursive function that should draw a node
    // We use a std::function because a simple lambda cannot be recursive
    const std::function<void(int, const glm::mat4 &)> drawNode =
        [&](int nodeIdx, const glm::mat4 &parentMatrix) {
          // Get the node and compute its modelMatrix.
          auto &node = model.nodes[nodeIdx];
          const glm::mat4 modelMatrix =
              getLocalToWorldMatrix(node, parentMatrix);

          if (node.mesh >= 0) {
            // Compute useful matrices.
            const glm::mat4 modelViewMatrix = viewMatrix * modelMatrix;
            const glm::mat4 modelViewProjectionMatrix =
                projMatrix * modelViewMatrix;
            const glm::mat4 normalMatrix =
                glm::transpose(glm::inverse(modelViewMatrix));

            // Send them to the GPU.
            glUniformMatrix4fv(modelViewMatrixLocation, 1, GL_FALSE,
                glm::value_ptr(modelViewMatrix));
            glUniformMatrix4fv(modelViewProjMatrixLocation, 1, GL_FALSE,
                glm::value_ptr(modelViewProjectionMatrix));
            glUniformMatrix4fv(normalMatrixLocation, 1, GL_FALSE,
                glm::value_ptr(normalMatrix));

            // Get the mesh and its corresponding vaoRange.
            auto &mesh = model.meshes[node.mesh];
            auto &vaoRange = meshIndexToVaoRange[node.mesh];

            // Draw each primitive of the mesh.
            for (auto primIdx = 0; primIdx < mesh.primitives.size();
                 ++primIdx) {
              auto &primitive = mesh.primitives[primIdx];
              auto vao = vertexArrayObjects[vaoRange.begin + primIdx];

              bindMaterial(primitive.material);

              glBindVertexArray(vao);

              // If the primitive uses IBO.
              if (primitive.indices >= 0) {
                auto &accessor = model.accessors[primitive.indices];
                auto &bufferView = model.bufferViews[accessor.bufferView];
                const auto byteOffset =
                    accessor.byteOffset + bufferView.byteOffset;

                glDrawElements(primitive.mode, GLsizei(accessor.count),
                    accessor.componentType, (const GLvoid *)byteOffset);
              } else { // If not using IBO.
                const auto accessorIdx = (*begin(primitive.attributes)).second;
                const auto &accessor = model.accessors[accessorIdx];

                glDrawArrays(primitive.mode, 0, GLsizei(accessor.count));
              }
            }
            // Unbind the VAO.
            glBindVertexArray(0);
          }

          // Recursively call the function on all children.
          for (const auto childIdx : node.children) {
            drawNode(childIdx, modelMatrix);
          }
        };

    // Draw the scene referenced by gltf file
    if (model.defaultScene >= 0) {
      // Draw all nodes.
      const auto &nodes = model.scenes[model.defaultScene].nodes;
      for (const auto &nodeIdx : nodes) {
        drawNode(nodeIdx, glm::mat4(1));
      }
    }
  };

  // Render the image in a PNG file.
  if (!m_OutputPath.empty()) {
    std::vector<unsigned char> pixels(m_nWindowWidth * m_nWindowHeight * 3);
    renderToImage(m_nWindowWidth, m_nWindowHeight, 3, pixels.data(),
        [&]() { drawScene(cameraController->getCamera()); });

    flipImageYAxis(m_nWindowWidth, m_nWindowHeight, 3, pixels.data());

    const auto strPath = m_OutputPath.string();
    stbi_write_png(
        strPath.c_str(), m_nWindowWidth, m_nWindowHeight, 3, pixels.data(), 0);

    return 0;
  }

  // Loop until the user closes the window
  for (auto iterationCount = 0u; !m_GLFWHandle.shouldClose();
       ++iterationCount) {
    const auto seconds = glfwGetTime();

    const auto camera = cameraController->getCamera();
    drawScene(camera);

    // GUI code:
    imguiNewFrame();

    {
      ImGui::Begin("GUI");
      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
          1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
      if (ImGui::CollapsingHeader("Camera", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Text("eye: %.3f %.3f %.3f", camera.eye().x, camera.eye().y,
            camera.eye().z);
        ImGui::Text("center: %.3f %.3f %.3f", camera.center().x,
            camera.center().y, camera.center().z);
        ImGui::Text(
            "up: %.3f %.3f %.3f", camera.up().x, camera.up().y, camera.up().z);

        ImGui::Text("front: %.3f %.3f %.3f", camera.front().x, camera.front().y,
            camera.front().z);
        ImGui::Text("left: %.3f %.3f %.3f", camera.left().x, camera.left().y,
            camera.left().z);

        if (ImGui::Button("CLI camera args to clipboard")) {
          std::stringstream ss;
          ss << "--lookat " << camera.eye().x << "," << camera.eye().y << ","
             << camera.eye().z << "," << camera.center().x << ","
             << camera.center().y << "," << camera.center().z << ","
             << camera.up().x << "," << camera.up().y << "," << camera.up().z;
          const auto str = ss.str();
          glfwSetClipboardString(m_GLFWHandle.window(), str.c_str());
        }

        // Deal with camera type.
        static int cameraType = 0;
        const auto trackballRadio =
            ImGui::RadioButton("Trackball", &cameraType, 0);
        ImGui::SameLine();
        const auto fpsRadio = ImGui::RadioButton("FPS", &cameraType, 1);

        const bool cameraTypeChanged = trackballRadio || fpsRadio;
        if (cameraTypeChanged) {
          const auto cam = cameraController->getCamera();
          if (cameraType == 0) {
            cameraController = std::make_unique<TrackballCameraController>(
                m_GLFWHandle.window(), 0.6f * maxDistance);
          } else {
            cameraController = std::make_unique<FirstPersonCameraController>(
                m_GLFWHandle.window(), 0.6f * maxDistance);
          }
          cameraController->setCamera(cam);
        }
      }
      if (ImGui::CollapsingHeader("Light", ImGuiTreeNodeFlags_DefaultOpen)) {
        // Add sliders for direction.
        static bool lightFromCamera = false;
        const auto lightFromCamCheckbox =
            ImGui::Checkbox("Light from Camera", &lightFromCamera);

        if (lightFromCamera) {
          lightDirection = glm::vec3(glm::inverse(camera.getViewMatrix()) *
                                     glm::vec4(0.f, 0.f, 1.f, 0.f));
        } else { // Light direction not computed from camera.
          static float theta = 0.f;
          static float phi = 0.f;

          // If the box was just unchecked, compute direction.
          if (lightFromCamCheckbox)
            lightDirection = glm::vec3(glm::sin(theta) * glm::cos(phi),
                glm::cos(theta), glm::sin(theta) * glm::sin(phi));

          const auto thetaSlider = ImGui::SliderFloat(
              "Light direction: theta", &theta, 0.f, 3.1415f, "%.2f");
          const auto phiSlider = ImGui::SliderFloat(
              "Light direction: phi", &phi, 0.f, 2 * 3.1415f, "%.2f");
          if (thetaSlider || phiSlider) {
            lightDirection = glm::vec3(glm::sin(theta) * glm::cos(phi),
                glm::cos(theta), glm::sin(theta) * glm::sin(phi));
          }
        }

        // Add color picker for intensity.
        static float color[3] = {1.f, 1.f, 1.f};
        static float intensity = 1.f;
        const auto colorPicker = ImGui::ColorEdit3("Light color", color);
        const auto intensityInput = ImGui::SliderFloat(
            "Light Intensity", &intensity, 0.0f, 10.0f, "%.2f");
        if (colorPicker || intensityInput) {
          lightIntensity = glm::vec3(color[0], color[1], color[2]) * intensity;
        }

        // Ambiant occlusion.
        const auto ambiantOcclusionCheckbox =
            ImGui::Checkbox("Ambiant Occlusion", &ambiantOcclusion);
      }
      ImGui::End();
    }

    imguiRenderFrame();

    glfwPollEvents(); // Poll for and process events

    auto ellapsedTime = glfwGetTime() - seconds;
    auto guiHasFocus =
        ImGui::GetIO().WantCaptureMouse || ImGui::GetIO().WantCaptureKeyboard;
    if (!guiHasFocus) {
      cameraController->update(float(ellapsedTime));
    }

    m_GLFWHandle.swapBuffers(); // Swap front and back buffers
  }

  // TODO clean up allocated GL data

  return 0;
}

ViewerApplication::ViewerApplication(const fs::path &appPath, uint32_t width,
    uint32_t height, const fs::path &gltfFile,
    const std::vector<float> &lookatArgs, const std::string &vertexShader,
    const std::string &fragmentShader, const fs::path &output) :
    m_nWindowWidth(width),
    m_nWindowHeight(height),
    m_AppPath{appPath},
    m_AppName{m_AppPath.stem().string()},
    m_ImGuiIniFilename{m_AppName + ".imgui.ini"},
    m_ShadersRootPath{m_AppPath.parent_path() / "shaders"},
    m_gltfFilePath{gltfFile},
    m_OutputPath{output}
{
  if (!lookatArgs.empty()) {
    m_hasUserCamera = true;
    m_userCamera =
        Camera{glm::vec3(lookatArgs[0], lookatArgs[1], lookatArgs[2]),
            glm::vec3(lookatArgs[3], lookatArgs[4], lookatArgs[5]),
            glm::vec3(lookatArgs[6], lookatArgs[7], lookatArgs[8])};
  }

  if (!vertexShader.empty()) {
    m_vertexShader = vertexShader;
  }

  if (!fragmentShader.empty()) {
    m_fragmentShader = fragmentShader;
  }

  ImGui::GetIO().IniFilename =
      m_ImGuiIniFilename.c_str(); // At exit, ImGUI will store its windows
                                  // positions in this file

  glfwSetKeyCallback(m_GLFWHandle.window(), keyCallback);

  printGLVersion();
}

bool ViewerApplication::loadGltfFile(tinygltf::Model &model) const
{
  tinygltf::TinyGLTF loader;
  std::string err;
  std::string warn;

  bool ret =
      loader.LoadASCIIFromFile(&model, &err, &warn, m_gltfFilePath.string());

  if (!warn.empty()) {
    std::cerr << "Warn: " << warn.c_str() << std::endl;
  }

  if (!err.empty()) {
    std::cerr << "Err: " << err.c_str() << std::endl;
  }

  return ret;
}

std::vector<GLuint> ViewerApplication::createBufferObjects(
    const tinygltf::Model &model) const
{
  // Initialize the identifiers vector.
  std::vector<GLuint> bufferObjects(model.buffers.size(), 0);
  // Generate the identifiers.
  glGenBuffers(GLsizei(bufferObjects.size()), bufferObjects.data());
  // Bind the data of each buffer to one identifier.
  for (size_t i = 0; i < bufferObjects.size(); ++i) {
    glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[i]);
    glBufferStorage(GL_ARRAY_BUFFER, model.buffers[i].data.size(),
        model.buffers[i].data.data(), 0);
  }
  // Unbind the GL_ARRAY_BUFFER.
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  return bufferObjects;
}

std::vector<GLuint> ViewerApplication::createVertexArrayObjects(
    const tinygltf::Model &model, const std::vector<GLuint> &bufferObjects,
    std::vector<VaoRange> &meshIndexToVaoRange)
{
  // Empty vector of VAO identifiers.
  std::vector<GLuint> vertexArrayObjects;

  // Loop though all the meshes of the model.
  for (size_t meshIdx = 0; meshIdx < model.meshes.size(); ++meshIdx) {
    // Get essential information.
    const auto &mesh = model.meshes[meshIdx];
    const auto nbPrimitives = mesh.primitives.size();
    const auto vaoOffset = vertexArrayObjects.size();

    // Create place for all the primitive's VAO identifiers of this mesh.
    vertexArrayObjects.resize(vaoOffset + nbPrimitives);
    // Store the offsets as it will be used during rendering.
    meshIndexToVaoRange.push_back(
        VaoRange{GLsizei(vaoOffset), GLsizei(nbPrimitives)});

    // Generate identifiers for all of those primitives.
    glGenVertexArrays(GLsizei(nbPrimitives), &vertexArrayObjects[vaoOffset]);

    // Loop through all the primitives of the mesh.
    for (size_t primIdx = 0; primIdx < nbPrimitives; ++primIdx) {
      // Get essential information.
      const auto &primitive = mesh.primitives[primIdx];
      const auto currentVao = vertexArrayObjects[vaoOffset + primIdx];

      // Bind this VAO.
      glBindVertexArray(currentVao);

      float *positions = nullptr;
      size_t positionsCount;

      float *texCoords = nullptr;
      size_t texCoordsCount;

      // Set all the attributes. Can be done better.
      { // POSITION
        const auto iterator = primitive.attributes.find("POSITION");
        if (iterator != end(primitive.attributes)) {
          const auto accessorIdx = (*iterator).second;
          const auto &accessor = model.accessors[accessorIdx];
          const auto &bufferView = model.bufferViews[accessor.bufferView];
          const auto bufferIdx = bufferView.buffer;

          // Enable the vertex attrib array corresponding to POSITION.
          glEnableVertexAttribArray(VERTEX_ATTRIB_POSITION_IDX);
          glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[bufferIdx]);

          const auto byteOffset = accessor.byteOffset + bufferView.byteOffset;
          glVertexAttribPointer(VERTEX_ATTRIB_POSITION_IDX, accessor.type,
              accessor.componentType, GL_FALSE, GLsizei(bufferView.byteStride),
              (const GLvoid *)byteOffset);

          const float *interpretPositions = reinterpret_cast<const float *>(
              &model.buffers[bufferIdx].data[byteOffset]);

          positions = (float *)(interpretPositions);
          positionsCount = accessor.count;
        }
      }
      { // NORMAL
        const auto iterator = primitive.attributes.find("NORMAL");
        if (iterator != end(primitive.attributes)) {
          const auto accessorIdx = (*iterator).second;
          const auto &accessor = model.accessors[accessorIdx];
          const auto &bufferView = model.bufferViews[accessor.bufferView];
          const auto bufferIdx = bufferView.buffer;

          // Enable the vertex attrib array corresponding to NORMAL.
          glEnableVertexAttribArray(VERTEX_ATTRIB_NORMAL_IDX);
          glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[bufferIdx]);

          const auto byteOffset = accessor.byteOffset + bufferView.byteOffset;
          glVertexAttribPointer(VERTEX_ATTRIB_NORMAL_IDX, accessor.type,
              accessor.componentType, GL_FALSE, GLsizei(bufferView.byteStride),
              (const GLvoid *)byteOffset);
        }
      }
      { // TEXCOORD_0
        const auto iterator = primitive.attributes.find("TEXCOORD_0");
        if (iterator != end(primitive.attributes)) {
          const auto accessorIdx = (*iterator).second;
          const auto &accessor = model.accessors[accessorIdx];
          const auto &bufferView = model.bufferViews[accessor.bufferView];
          const auto bufferIdx = bufferView.buffer;

          // Enable the vertex attrib array corresponding to TEXCOORD0.
          glEnableVertexAttribArray(VERTEX_ATTRIB_TEXCOORD0_IDX);
          glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[bufferIdx]);

          const auto byteOffset = accessor.byteOffset + bufferView.byteOffset;
          glVertexAttribPointer(VERTEX_ATTRIB_TEXCOORD0_IDX, accessor.type,
              accessor.componentType, GL_FALSE, GLsizei(bufferView.byteStride),
              (const GLvoid *)byteOffset);

          const float *interpretTexCoords = reinterpret_cast<const float *>(
              &model.buffers[bufferIdx].data[byteOffset]);

          texCoords = (float *)(interpretTexCoords);
          texCoordsCount = accessor.count;
        }
      }

      // Set the index buffer if there is one.
      if (primitive.indices >= 0) {
        const auto &accessor = model.accessors[primitive.indices];
        const auto &bufferView = model.bufferViews[accessor.bufferView];
        const auto bufferIdx = bufferView.buffer;
        const auto byteOffset = accessor.byteOffset + bufferView.byteOffset;

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferObjects[bufferIdx]);

        // Compute tangents for each vertex using the indices.
        std::vector<glm::vec3> tangents(positionsCount);
        const auto *indices = reinterpret_cast<const unsigned short *>(
            &model.buffers[bufferIdx].data[byteOffset]);

        // Algo from http://ogldev.atspace.co.uk/www/tutorial26/tutorial26.html
        for (auto i = 0; i < accessor.count; i += 3) {
          // Vertex indices.
          const auto vertex0 = indices[i];
          const auto vertex1 = indices[i + 1];
          const auto vertex2 = indices[i + 2];

          // Vertex 0.
          auto indexPos = vertex0 * 3;
          auto indexTex = vertex0 * 2;
          const auto posV0 = glm::vec3(positions[indexPos + 0],
              positions[indexPos + 1], positions[indexPos + 2]);
          const auto texV0 =
              glm::vec2(texCoords[indexTex + 0], texCoords[indexTex + 1]);
          // Vertex 1.
          indexPos = vertex1 * 3;
          indexTex = vertex1 * 2;
          const auto posV1 = glm::vec3(positions[indexPos + 0],
              positions[indexPos + 1], positions[indexPos + 2]);
          const auto texV1 =
              glm::vec2(texCoords[indexTex + 0], texCoords[indexTex + 1]);
          // Vertex 2.
          indexPos = vertex2 * 3;
          indexTex = vertex2 * 2;
          const auto posV2 = glm::vec3(positions[indexPos + 0],
              positions[indexPos + 1], positions[indexPos + 2]);
          const auto texV2 =
              glm::vec2(texCoords[indexTex + 0], texCoords[indexTex + 1]);

          // Compute edges.
          glm::vec3 edge1 = posV1 - posV0;
          glm::vec3 edge2 = posV2 - posV0;

          // Compute deltas.
          float deltaU1 = texV1.x - texV0.x;
          float deltaV1 = texV1.y - texV0.y;
          float deltaU2 = texV2.x - texV0.x;
          float deltaV2 = texV2.y - texV0.y;

          float f = 1.0f / (deltaU1 * deltaV2 - deltaU2 * deltaV1);

          glm::vec3 tangent;
          tangent.x = f * (deltaV2 * edge1.x - deltaV1 * edge2.x);
          tangent.y = f * (deltaV2 * edge1.y - deltaV1 * edge2.y);
          tangent.z = f * (deltaV2 * edge1.z - deltaV1 * edge2.z);

          // Add the tangent to vertices tangent.
          tangents[vertex0] += tangent;
          tangents[vertex1] += tangent;
          tangents[vertex2] += tangent;
        }

        // Normalize every tangent.
        for (auto &t : tangents) {
          t = glm::normalize(t);
        }

        GLuint bufObj;
        glGenBuffers(1, &bufObj);
        glBindBuffer(GL_ARRAY_BUFFER, bufObj);
        glBufferStorage(GL_ARRAY_BUFFER, tangents.size() * sizeof(glm::vec3),
            tangents.data(), 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        tangentBuffers.push_back(bufObj); // Used only to clean up the buffers.

        // Enable the vertex attrib array corresponding to TANGENT.
        glEnableVertexAttribArray(VERTEX_ATTRIB_TANGENT_IDX);
        glBindBuffer(GL_ARRAY_BUFFER, bufObj);
        glVertexAttribPointer(VERTEX_ATTRIB_TANGENT_IDX, 3, GL_FLOAT, GL_FALSE,
            sizeof(float) * 3, 0);
      }
    }
  }

  glBindVertexArray(0);
  return vertexArrayObjects;
}

std::vector<GLuint> ViewerApplication::createTextureObjects(
    const tinygltf::Model &model) const
{
  // Texture identifiers.
  std::vector<GLuint> textureIds(model.textures.size(), 0);
  glGenTextures(GLsizei(textureIds.size()), textureIds.data());

  // Default sampler.
  tinygltf::Sampler defaultSampler;
  defaultSampler.minFilter = GL_LINEAR;
  defaultSampler.magFilter = GL_LINEAR;
  defaultSampler.wrapS = GL_REPEAT;
  defaultSampler.wrapT = GL_REPEAT;
  defaultSampler.wrapR = GL_REPEAT;

  // Loop through each texture.
  for (auto i = 0; i < textureIds.size(); ++i) {
    // Get texture.
    const auto &texture = model.textures[i];
    // Get the corresponding image. Make sure there is one.
    assert(texture.source >= 0);
    const auto &image = model.images[texture.source];
    // Get the corresponding sampler (or the default one if there is not).
    const auto &sampler =
        texture.sampler >= 0 ? model.samplers[texture.sampler] : defaultSampler;

    // Bind texture object to target GL_TEXTURE_2D.
    glBindTexture(GL_TEXTURE_2D, textureIds[i]);
    // Set image data.
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width, image.height, 0,
        GL_RGBA, image.pixel_type, image.image.data());
    // Set sampling parameters.
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        sampler.minFilter != -1 ? sampler.minFilter : GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
        sampler.magFilter != -1 ? sampler.magFilter : GL_LINEAR);
    // Set wrapping parameters.
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, sampler.wrapS);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, sampler.wrapT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, sampler.wrapR);

    // Check if mipmapping is needed.
    if (sampler.minFilter == GL_NEAREST_MIPMAP_NEAREST ||
        sampler.minFilter == GL_NEAREST_MIPMAP_LINEAR ||
        sampler.minFilter == GL_LINEAR_MIPMAP_NEAREST ||
        sampler.minFilter == GL_LINEAR_MIPMAP_LINEAR) {
      glGenerateMipmap(GL_TEXTURE_2D);
    }
  }
  // Unbind texture.
  glBindTexture(GL_TEXTURE_2D, 0);

  return textureIds;
}