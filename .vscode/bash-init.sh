clone_gltf_samples() {
  git clone https://github.com/KhronosGroup/glTF-Sample-Models gltf-sample-models
}

get_sponza_sample() {
  mkdir -p gltf-sample-models/2.0/ || true
  curl -L https://gltf-viewer-tutorial.gitlab.io/assets/Sponza.zip > gltf-sample-models/2.0/Sponza.zip
  unzip gltf-sample-models/2.0/Sponza.zip -d gltf-sample-models/2.0/
  rm gltf-sample-models/2.0/Sponza.zip
}

get_damaged_helmet_sample() {
  mkdir -p gltf-sample-models/2.0/ || true
  curl -L https://gltf-viewer-tutorial.gitlab.io/assets/DamagedHelmet.zip > gltf-sample-models/2.0/DamagedHelmet.zip
  unzip gltf-sample-models/2.0/DamagedHelmet.zip -d gltf-sample-models/2.0/
  rm gltf-sample-models/2.0/DamagedHelmet.zip
}

cmake_clean() {
  rm -rf build
  rm -rf dist
}

cmake_prepare() {
  # Note: $@ forwards argument so you can add additional arguments
  # e.g: cmake_prepare -DCMAKE_BUILD_TYPE=Release to configure the build solution in release mode with gcc
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX=./dist $@
}

cmake_build() {
  # Note: $@ forwards argument so you can add additional arguments
  # e.g: cmake_build --config release to build in release mode with Visual Studio Compiler
  cmake --build build -j $@
}

cmake_install() {
  # Note: $@ forwards argument so you can add additional arguments
  # e.g: cmake_install --config release to build in release mode, then install, with Visual Studio Compiler
  cmake --build build -j --target install $@
}

view_sponza() {
  cmake_prepare
  cmake_install
  dist/gltf-viewer viewer gltf-sample-models/2.0/Sponza/glTF/Sponza.gltf \
    --lookat -5.26056,6.59932,0.85661,-4.40144,6.23486,0.497347,0.342113,0.931131,-0.126476
}

view_helmet() {
  cmake_prepare
  cmake_install
  dist/gltf-viewer viewer gltf-sample-models/2.0/DamagedHelmet/glTF/DamagedHelmet.gltf
}

render_sponza() {
  cmake_prepare
  cmake_install
  [ ! -d output-images ] && mkdir output-images
  dist/gltf-viewer viewer gltf-sample-models/2.0/Sponza/glTF/Sponza.gltf \
    --lookat -5.26056,6.59932,0.85661,-4.40144,6.23486,0.497347,0.342113,0.931131,-0.126476 \
    --output output-images/sponza.png
}

render_helmet() {
  cmake_prepare
  cmake_install
  dist/gltf-viewer viewer gltf-sample-models/2.0/DamagedHelmet/glTF/DamagedHelmet.gltf \
    --output output-images/helmet.png
}
