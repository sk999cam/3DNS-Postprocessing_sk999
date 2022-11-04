echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v212\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v212\fluent\ntbin\win64\tell.exe" WL-CH-009 63220 CLEANUP_EXITING
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 3884) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 14728) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 23616)
del "E:\MATLAB\RANS_optimisation\cleanup-fluent-WL-CH-009-14728.bat"
