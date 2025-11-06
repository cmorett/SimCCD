# SimCCD_Log_Complete 

Dependencias
- Geant4 (probado con geant4-11-03-patch-02)
- Assimp (importador COLLADA) — instala vía vcpkg
- CMake (probaado con 4.1.2)
- Visual Studio / MSVC (x64) en Windows
- vcpkg

Créditos
- Assimp contributors — importador COLLADA.

Enlace al archivo Onshape (tab 1309 option Copy 1)
- Enlace Onshape: https://cad.onshape.com/documents/9783a63d985190bd76485a0b/w/865405b6c01c66e03051095c/e/07f8aa569993da8d219a7d36?renderMode=0&uiState=690beb6c005da20a7e5f1127

CMake: variable requerida (vcpkg)
- Antes de configurar, pasar el toolchain de vcpkg a CMake:
  -DCMAKE_TOOLCHAIN_FILE="C:\vcpkg\scripts\buildsystems\vcpkg.cmake"

Ejemplo de comandos (Windows, Release)
- cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE="C:\vcpkg\scripts\buildsystems\vcpkg.cmake" -DCMAKE_BUILD_TYPE=Release
- cmake --build build --config Release

Nota
- Usar la configuración Release (coincidente con las librerías Geant4 precompiladas) para evitar errores en Visual Studio (para Windows).
- Preferible mantener la información de unidades correcta en el .dae (meter="1") y convertir en el loader.
- Para cambios geométricos finos se puede reexportar/normalizar en Blender (importar, escalar, aplicar transformaciones y exportar).

