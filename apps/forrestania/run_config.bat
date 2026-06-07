@echo off
setlocal
set EXE=..\..\build\release\forrestania_invert.exe
set CONFIG_DIR=%~dp0

if "%1"=="" (
    echo Usage: run_config.bat ^<name^>
    echo.
    echo Available configs:
    echo   standard    - 27 split meshes from GMM clustering ^(v2 cleaned^)
    echo   baseline    -  5 gated clusters, baseline properties
    echo   perturbed   -  5 gated clusters, denser + more magnetic target
    echo.
    echo Examples:
    echo   run_config standard
    echo   run_config perturbed
    exit /b 1
)

if /i "%1"=="standard" (
    set CONFIG=%CONFIG_DIR%resolved_config.ini
    echo === STANDARD: 27 split meshes, GMM clustering ===
)
if /i "%1"=="baseline" (
    set CONFIG=%CONFIG_DIR%User_Perturbed_Physical_Properties\baseline\resolved_config.ini
    echo === BASELINE: 5 gated clusters, baseline properties ===
)
if /i "%1"=="perturbed" (
    set CONFIG=%CONFIG_DIR%User_Perturbed_Physical_Properties\perturbed\resolved_config.ini
    echo === PERTURBED: 5 gated clusters, denser target ===
)

if not defined CONFIG (
    echo Unknown config: %1
    exit /b 1
)

echo Config: %CONFIG%
echo.
"%EXE%" "%CONFIG%"

