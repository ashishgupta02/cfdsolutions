# CFD Solutions

A suite of computational fluid dynamics libraries and solvers (2D/3D Euler and
multi-phase flow), mesh generation/quality tools, and supporting numerical
libraries. The project is built with **CMake**.

## Repository layout

All build files live under `Trunk/`:

| Component | Type | Notes |
|-----------|------|-------|
| `UTILS` | static lib | core utilities |
| `MATH` | static lib | linear algebra, geometry |
| `CFDDB` | static lib | CFD database / mesh I/O (optional NetCDF) |
| `NISTThermo` | static lib | NIST thermodynamics (Fortran) |
| `EOS` | static lib | equation-of-state models |
| `VERDICT` | static lib | mesh quality metrics |
| `mpflow3d` | executable | multi-phase 3D flow solver |
| `euler2d` / `euler3d` | executable | Euler solvers |
| `mesh2d` | executable | 2D mesh generator |
| `ugrid2su2` | executable | grid converter |
| `par_euler3d` | executable | parallel Euler solver (**requires MPI**) |
| `meshqa` | executable | mesh QA (**requires CGNS**) |
| `postcfd` | executable | post-processing (**requires CGNS**) |

## Prerequisites

- **CMake ≥ 3.16**
- **Ninja** (or Unix Makefiles). A Fortran-capable generator is required because
  `NISTThermo` is written in Fortran — **the Visual Studio generator does not
  support Fortran and cannot be used.**
- A supported compiler toolchain:
  - **GNU**: `gcc`, `g++`, `gfortran`
  - **Intel oneAPI**: `icx`, `icpx`, `ifx`
  - **Clang**: `clang`, `clang++` (Fortran sources are compiled with the
    auto-detected Fortran compiler, typically `gfortran`)
- Optional external dependencies, enabled by the options below:
  - **MPI** (e.g. MPICH/Open MPI/Intel MPI) — for `par_euler3d` and parallel solver paths
  - **NetCDF** — for the CFDDB NetCDF I/O path
  - **CGNS** — for `meshqa` and `postcfd`

## Build options

| Option | Default | Effect |
|--------|---------|--------|
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` or `Release` |
| `ENABLE_MPI` | `OFF` | Build MPI paths and the `par_euler3d` target |
| `ENABLE_NETCDF` | `OFF` | Compile the CFDDB NetCDF I/O source (`netcdfIO.c`) |
| `ENABLE_CGNS` | `OFF` | Build the CGNS-dependent `meshqa` and `postcfd` targets |
| `ENABLE_NATIVE_ARCH` | `ON` | Add `-march=native` (tunes for the build machine's CPU) |
| `ENABLE_STRICT_WARNINGS` | `ON` | `-Wall -Wextra -Wpedantic` on GNU/Clang |
| `BUILD_SHARED_LIBS` | `OFF` | Build the internal libraries as shared objects |

## Quick start (CMake presets)

Presets are defined in `Trunk/CMakePresets.json` for every supported combination
of **compiler × build type × MPI**. All presets use the Ninja generator and place
their build tree in `Trunk/build/<preset>`.

```bash
cd Trunk

# Configure and build a GNU Release build
cmake --preset gnu-release
cmake --build --preset gnu-release
```

Available presets:

| Compiler | Debug | Release | Debug + MPI | Release + MPI |
|----------|-------|---------|-------------|---------------|
| GNU | `gnu-debug` | `gnu-release` | `gnu-debug-mpi` | `gnu-release-mpi` |
| Intel | `intel-debug` | `intel-release` | `intel-debug-mpi` | `intel-release-mpi` |
| Clang | `clang-debug` | `clang-release` | `clang-debug-mpi` | `clang-release-mpi` |

List them at any time with:

```bash
cd Trunk
cmake --list-presets
```

## Manual configuration (without presets)

If you prefer to drive CMake directly, always pass a Fortran-capable generator
(`-G Ninja`) and select the compilers explicitly.

### GNU

```bash
cd Trunk
cmake -S . -B build/gnu -G Ninja \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/gnu
```

### Intel oneAPI

```bash
# Source the Intel environment first (setvars.sh / setvars.bat)
cd Trunk
cmake -S . -B build/intel -G Ninja \
  -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/intel
```

### Clang

```bash
cd Trunk
cmake -S . -B build/clang -G Ninja \
  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_BUILD_TYPE=Debug
cmake --build build/clang
```

### Enabling optional features

```bash
# MPI (also builds par_euler3d)
cmake -S . -B build/mpi -G Ninja -DENABLE_MPI=ON

# CGNS tools (meshqa, postcfd)
cmake -S . -B build/cgns -G Ninja -DENABLE_CGNS=ON

# NetCDF I/O in CFDDB
cmake -S . -B build/netcdf -G Ninja -DENABLE_NETCDF=ON

# Portable binaries (disable -march=native)
cmake -S . -B build/portable -G Ninja -DENABLE_NATIVE_ARCH=OFF
```

## Windows notes

On Windows the recommended setup is a GNU toolchain (e.g. MSYS2 UCRT64 providing
`gcc`/`g++`/`gfortran`) together with Ninja. Make sure `cmake`, `ninja`, and the
compilers are on `PATH`, then use the presets exactly as above. The Visual Studio
generator will not work because of the Fortran component.

## Compiler flags

All optimization, floating-point, and warning flags are defined centrally in
`Trunk/cmake/CompilerFlags.cmake` (there are no per-target flag overrides). The
Release profile mirrors the project's historical tuning:
`-O3 -march=native -ffast-math -funroll-all-loops -finline-limit=150000` and the
related fast-math switches for GNU (with equivalent mappings for Clang and Intel).

> **Note:** `-ffast-math` implies `-ffinite-math-only`, which lets the compiler
> assume values are never `NaN`/`Inf`. `isnan()`-based divergence checks may be
> optimized away in Release builds. Build `Debug` (or disable the fast-math flags)
> when you need reliable NaN detection.
