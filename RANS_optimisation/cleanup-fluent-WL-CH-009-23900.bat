echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v212\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v212\fluent\ntbin\win64\tell.exe" WL-CH-009 54392 CLEANUP_EXITING
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 29756) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 23900) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 19796)
del "E:\MATLAB\RANS_optimisation\cleanup-fluent-WL-CH-009-23900.bat"
