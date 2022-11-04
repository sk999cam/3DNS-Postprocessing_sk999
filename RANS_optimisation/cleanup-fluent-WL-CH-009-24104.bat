echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v212\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v212\fluent\ntbin\win64\tell.exe" WL-CH-009 62903 CLEANUP_EXITING
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 25512) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 24104) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 26068)
del "E:\MATLAB\RANS_optimisation\cleanup-fluent-WL-CH-009-24104.bat"
