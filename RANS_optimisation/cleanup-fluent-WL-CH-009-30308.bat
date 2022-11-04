echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v212\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v212\fluent\ntbin\win64\tell.exe" WL-CH-009 63060 CLEANUP_EXITING
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 7124) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 30308) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 29968)
del "E:\MATLAB\RANS_optimisation\cleanup-fluent-WL-CH-009-30308.bat"
