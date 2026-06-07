@echo off
REM LithoInvert3D — Full Build Script
REM Requires: Visual Studio 2022, Qt 6.2+ with qmake, Eigen 3.4.0 in vendor/eigen/
REM
REM Build order (dependency chain):
REM   litho-core → litho-surface → litho-model → litho-forward/litho-em/litho-io/litho-regularization → litho-inversion

setlocal enabledelayedexpansion

REM Locate vcvars64.bat
if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
) else if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" (
    call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat"
) else (
    echo ERROR: vcvars64.bat not found. Run from Developer Command Prompt for VS 2022.
    exit /b 1
)

set "ROOT=%~dp0.."
set "BUILD_DIR=%ROOT%\build\release"
set "EIGEN_DIR=%ROOT%\vendor\eigen"

if not exist "%EIGEN_DIR%\Eigen" (
    echo ERROR: Eigen not found at %EIGEN_DIR%
    echo Copy Eigen 3.4.0 to vendor/eigen/
    exit /b 1
)

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"

REM Module build order (respects dependency chain)
set MODULES=litho-core litho-surface litho-model litho-forward litho-em litho-io litho-regularization litho-inversion

for %%M in (%MODULES%) do (
    echo.
    echo ========================================
    echo Building %%M...
    echo ========================================
    cd "%ROOT%\modules\%%M"
    qmake "%%M.pro"
    if errorlevel 1 (
        echo ERROR: qmake failed for %%M
        exit /b 1
    )
    nmake release
    if errorlevel 1 (
        echo ERROR: nmake failed for %%M
        exit /b 1
    )
    REM Copy .lib to build/release/ for downstream linking
    if exist "release\%%M.lib" (
        copy /Y "release\%%M.lib" "%BUILD_DIR%\"
    )
)

echo.
echo ========================================
echo Building executable (apps/forrestania)...
echo ========================================
cd "%ROOT%\apps\forrestania"
qmake forrestania.pro
if errorlevel 1 (
    echo ERROR: qmake failed for forrestania
    exit /b 1
)
nmake release
if errorlevel 1 (
    echo ERROR: nmake failed for forrestania
    exit /b 1
)
copy /Y "release\forrestania_invert.exe" "%BUILD_DIR%\"

echo.
echo ========================================
echo BUILD COMPLETE
echo Executable: %BUILD_DIR%\forrestania_invert.exe
echo ========================================

