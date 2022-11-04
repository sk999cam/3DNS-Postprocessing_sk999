echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v212\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v212\fluent\ntbin\win64\tell.exe" WL-CH-009 62654 CLEANUP_EXITING
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 12604) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 26352) 
if /i "%LOCALHOST%"=="WL-CH-009" (%KILL_CMD% 12464)
del "E:\MATLAB\RANS_optimisation\cleanup-fluent-WL-CH-009-26352.bat"
