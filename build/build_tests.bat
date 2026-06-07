@echo off
setlocal
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" > nul
set "ROOT=%~dp0.."

echo ========================================
echo Building forrestania tests...
echo ========================================
cd "%ROOT%\apps\forrestania\tests"
qmake tests.pro
nmake release
if errorlevel 1 (
    echo ERROR: nmake failed for forrestania tests
    exit /b 1
)
echo Tests built successfully.
