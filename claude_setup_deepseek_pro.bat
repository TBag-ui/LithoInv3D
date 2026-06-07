@echo off
setlocal enabledelayedexpansion
title Claude Code Terminal Setup

:: ---- Store DeepSeek API key securely (one-time prompt) ----
set "KEYFILE=%USERPROFILE%\.deepseek_apikey"
if not exist "%KEYFILE%" (
    echo.
    echo DeepSeek API key not found. You only need to enter it once.
    set /p "DSKEY=Enter your DeepSeek API key: "
    echo !DSKEY! > "%KEYFILE%"
    echo Key saved for future use.
) else (
    set /p DSKEY=<"%KEYFILE%"
)

:: ---- Create a temporary batch file for the DeepSeek terminal ----
:: This avoids exposing the API key in the command-line of the new window
set "DSTMP=%TEMP%\claude_deepseek_env.bat"
(
echo @echo off
echo set ANTHROPIC_BASE_URL=https://api.deepseek.com/anthropic
echo set ANTHROPIC_API_KEY=%DSKEY%
echo set ANTHROPIC_MODEL=deepseek-v4-pro
echo set ANTHROPIC_DEFAULT_OPUS_MODEL=deepseek-v4-pro
echo set ANTHROPIC_DEFAULT_SONNET_MODEL=deepseek-v4-pro
echo set ANTHROPIC_DEFAULT_HAIKU_MODEL=deepseek-v4-flash
echo set CLAUDE_CODE_SUBAGENT_MODEL=deepseek-v4-flash
echo set CLAUDE_CODE_EFFORT_LEVEL=max
echo cd /d "%CD%"
echo echo.
echo echo ============================================
echo echo   DeepSeek PRO environment is active.
echo echo   Main tasks will use the powerful deepseek-v4-pro model.
echo echo   Sub-tasks will use the fast deepseek-v4-flash model.
echo echo   Run 'claude' to start.
echo echo ============================================
echo echo.
echo cmd /k
) > "%DSTMP%"

:: ---- Launch the DeepSeek terminal ----
start "Claude Code - DeepSeek" "%DSTMP%"

:: ---- Launch the Pro subscription terminal ----
start "Claude Code - Pro" cmd /k "cd /d %CD% && echo ============================================ && echo   Pro subscription terminal. && echo   Make sure you've run 'claude /login' once. && echo ============================================ && echo. && cmd /k"

:: ---- Cleanup (optional, the temp file is tiny) ----
:: del "%DSTMP%" 2>nul

echo.
echo Two terminals opened:
echo   [Claude Code - DeepSeek Pro] – ready with DeepSeek API (v4-pro main, flash sub-tasks)
echo   [Claude Code - Pro]      – ready for your Claude subscription
echo.
pause
