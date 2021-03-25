#pragma once

#include "cameras.hpp"
#include <glm/glm.hpp>

class CameraController
{
public:
  CameraController() :
      m_camera{glm::vec3(0), glm::vec3(0, 0, -1), glm::vec3(0, 1, 0)}
  {
  }
  virtual ~CameraController() = default;

  const Camera &getCamera() const { return m_camera; }
  void setCamera(const Camera &camera) { m_camera = camera; }

  // Must be implemented in subclasses.
  virtual bool update(float elapsedTime) = 0;

protected:
  // Current camera.
  Camera m_camera;
};

class FirstPersonCameraController : public CameraController
{
public:
  FirstPersonCameraController(GLFWwindow *window, float speed = 1.f,
      const glm::vec3 &worldUpAxis = glm::vec3(0, 1, 0)) :
      CameraController(),
      m_pWindow(window),
      m_fSpeed(speed),
      m_worldUpAxis(worldUpAxis)
  {
  }

  // Controller attributes, if put in a GUI, should be adapted
  void setSpeed(float speed) { m_fSpeed = speed; }

  float getSpeed() const { return m_fSpeed; }

  void increaseSpeed(float delta)
  {
    m_fSpeed += delta;
    m_fSpeed = glm::max(m_fSpeed, 0.f);
  }

  const glm::vec3 &getWorldUpAxis() const { return m_worldUpAxis; }

  void setWorldUpAxis(const glm::vec3 &worldUpAxis)
  {
    m_worldUpAxis = worldUpAxis;
  }

  // Update the view matrix based on input events and elapsed time
  // Return true if the view matrix has been modified
  bool update(float elapsedTime) override;

private:
  GLFWwindow *m_pWindow = nullptr;
  float m_fSpeed = 0.f;
  glm::vec3 m_worldUpAxis;

  // Input event state
  bool m_LeftButtonPressed = false;
  glm::dvec2 m_LastCursorPosition;
};

// todo Blender like camera
class TrackballCameraController : public CameraController
{
public:
  TrackballCameraController(GLFWwindow *window, float speed = 1.f,
      const glm::vec3 &worldUpAxis = glm::vec3(0, 1, 0)) :
      CameraController(),
      m_pWindow(window),
      m_fSpeed(speed),
      m_worldUpAxis(worldUpAxis)
  {
  }

  // Controller attributes, if put in a GUI, should be adapted
  void setSpeed(float speed) { m_fSpeed = speed; }

  float getSpeed() const { return m_fSpeed; }

  void increaseSpeed(float delta)
  {
    m_fSpeed += delta;
    m_fSpeed = glm::max(m_fSpeed, 0.f);
  }

  const glm::vec3 &getWorldUpAxis() const { return m_worldUpAxis; }

  void setWorldUpAxis(const glm::vec3 &worldUpAxis)
  {
    m_worldUpAxis = worldUpAxis;
  }

  // Update the view matrix based on input events and elapsed time
  // Return true if the view matrix has been modified
  bool update(float elapsedTime) override;

private:
  GLFWwindow *m_pWindow = nullptr;
  float m_fSpeed = 0.f;
  glm::vec3 m_worldUpAxis;

  // Input event state
  bool m_MiddleButtonPressed = false;
  glm::dvec2 m_LastCursorPosition;
};