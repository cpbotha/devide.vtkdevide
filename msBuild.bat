@echo off
REM build batch script for msdev under windows

REM this is how we did it with VC6
REM MSDEV.com VTKDEVIDE.dsw /MAKE "ALL_BUILD - Win32 RelWithDebInfo"

REM this is how we do it with VS.NET 2003 (VS71)
REM devenv ITK.sln /project ALL_BUILD /projectconfig "Debug|Win32" /build Debug
devenv ITK.sln /project ALL_BUILD /projectconfig "RelWithDebInfo|Win32" /build RelWithDebInfo




