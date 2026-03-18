@echo off
setlocal

REM ------------------------------------------------------------------------------
REM Author: Francisco Garre-Frutos
REM ------------------------------------------------------------------------------

set "ROOT=%~dp0.."
set "MANUSCRIPT_DIR=%ROOT%\output\manuscript"

if not exist "%MANUSCRIPT_DIR%" (
  mkdir "%MANUSCRIPT_DIR%"
)

for %%F in (
  "%ROOT%\*.aux"
  "%ROOT%\*.bbl"
  "%ROOT%\*.bcf"
  "%ROOT%\*.blg"
  "%ROOT%\*.fdb_latexmk"
  "%ROOT%\*.fls"
  "%ROOT%\*.log"
  "%ROOT%\*.out"
  "%ROOT%\*.run.xml"
  "%ROOT%\*.synctex.gz"
  "%ROOT%\*.tex"
  "%ROOT%\*.toc"
  "%ROOT%\*.xdv"
) do (
  if exist "%%~fF" (
    move /Y "%%~fF" "%MANUSCRIPT_DIR%" >nul
  )
)

endlocal
