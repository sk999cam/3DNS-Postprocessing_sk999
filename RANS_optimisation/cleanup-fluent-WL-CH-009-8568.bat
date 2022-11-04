echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v212\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v212\fluent\ntbin\win64\tell.exe" WL-CH-009 55440 CLEANUP_EXITING
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 14940) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 8568) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 20000)
del "E:\MATLAB\RANS_optimisation\cleanup-fluent-WL-CH-009-8568.bat"
