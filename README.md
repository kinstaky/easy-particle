# easy-particle

A C++17 library for basic nuclear-particle kinematics and energy-loss calculations.

## Features

- Particle model with relativistic energy/momentum updates.
- Two-body `Scatter` and two-fragment `Breakup` kinematics helpers.
- Material helpers built on `catima` (`GasMaterial`, `SolidMaterial`).
- Nuclear mass lookup from AMDC data file (`data/amdc_ion_2020.txt`) or optional hardcoded table.
- Path resolution that works in standalone builds and when consumed by another CMake project.

## Requirements

- CMake >= 3.10
- C++17 compiler
- cern ROOT

Dependencies fetched automatically by CMake:

- `catima`
- `cxxopts`

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

Build with tests:

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Data File Resolution

`PathManager::DataPath(filename)` searches candidate directories and returns the first directory that contains `filename`.

Search includes:

- executable-relative paths (`exec`, `..`, `../share/easy-particle`, ...)
- CMake-injected hints:
  - `EASY_PARTICLE_DATA_DIR`
  - `EASY_PARTICLE_SOURCE_DATA_DIR`

`nuclear_data.cpp` requests `"amdc_ion_2020.txt"` explicitly.

## Use In Another Project

### Option 1: Git submodule + `add_subdirectory`

Add repository as submodule (example path `third_party/easy-particle`) and in parent `CMakeLists.txt`:

```cmake
add_subdirectory(third_party/easy-particle)

target_link_libraries(your_target PRIVATE particle)
```

### Option 2: CMake `FetchContent`

```cmake
include(FetchContent)

FetchContent_Declare(
  easy_particle
  GIT_REPOSITORY <your-easy-particle-repo-url>
  GIT_TAG main
)
FetchContent_MakeAvailable(easy_particle)

target_link_libraries(your_target PRIVATE particle)
```

### Option 3: Installed package style (basic)

```bash
cmake -S . -B build
cmake --build build -j
cmake --install build --prefix /your/prefix
```

Then expose include/library paths from `/your/prefix` in your parent project.

## Minimal C++ Example

```cpp
#include "include/particle.h"

int main() {
  easyparticle::Particle beam(6, 16, 200.0);
  easyparticle::Particle target(1, 1);
  easyparticle::Particle f0(6, 16);
  easyparticle::Particle f1(1, 1);

  easyparticle::Scatter(
    beam,
    target,
    ROOT::Math::XYZVector(0.2, 0.0, 0.98).Unit(),
    f0,
    f1
  );

  return 0;
}
```

## Notes For Parent Projects

- The library is designed to be included as a subproject.
- Data path hints are provided via compile definitions from easy-particle's own CMake.
- If your packaging layout is custom, keep `amdc_ion_2020.txt` in one of the searched directories or patch candidate paths in `src/path_manager.cpp`.
