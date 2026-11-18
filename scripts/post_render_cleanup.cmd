@echo off
setlocal

set "ROOT=%~dp0.."
set "TARGET=%ROOT%\output\manuscript"

if not exist "%TARGET%" (
  mkdir "%TARGET%"
)

for %%B in (vmac_reversal_manuscript vmac_reversal_preprint vmac_reversal_supplementary) do (
  for %%E in (aux log tex xdv fdb_latexmk fls knit.md) do (
    if exist "%ROOT%\%%B.%%E" (
      move /Y "%ROOT%\%%B.%%E" "%TARGET%\" >nul
    )
  )
  if exist "%ROOT%\%%B_files" (
    move /Y "%ROOT%\%%B_files" "%TARGET%\" >nul
  )
)

endlocal
exit /b 0
