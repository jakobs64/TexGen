; Set the current version, this should change at each release
!define VERSION "3.14.0"
!define PYTHONDLL python312.dll
!define PYTHONVER "3.12"
!include "LogicLib.nsh"
!include "MUI.nsh"
!include "x64.nsh"

Var PythonDir

; The name of the installer
Name "TexGen Python3 ${VERSION}"

; The file to write
OutFile "texgen-Python312-${VERSION}.exe"

; Build a unicode installer
Unicode true

; TexGen version
VIProductVersion "${VERSION}"

; License file
; LicenseData "..\gpl.txt"

; The default installation directory
InstallDir $PROGRAMFILES64\TexGen

; Registry key to check for directory (so if you install again, it will
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\TexGen" "Install_Dir"

; Check that the python dll is present in the system path
Function CheckPython
; Set to 64-bit registry view
  SetRegView 64
  ClearErrors
  ReadRegStr $PythonDir HKLM SOFTWARE\Python\PythonCore\${PYTHONVER}\InstallPath ""
  ${If} ${Errors}
    ClearErrors
    ReadRegStr $PythonDir HKCU SOFTWARE\Python\PythonCore\${PYTHONVER}\InstallPath ""
  ${EndIf}

  ${If} ${Errors}
    MessageBox MB_YESNO|MB_ICONEXCLAMATION "Python ${PYTHONVER} doesn't seem to be installed, it is required for TexGen to run! \
                     $\nPlease download it from http://www.python.org/ \
                     $\n$\nContinue with installation anyway?" IDYES yes IDNO no
    no:
      Quit
    yes:
  ${EndIf}
  ; Set back to previous registry view
  SetRegView LastUsed
FunctionEnd

Function PythonPage
  Call CheckPython
  !insertmacro MUI_HEADER_TEXT "$(TEXT_PYTHON_TITLE)" "$(TEXT_PYTHON_SUBTITLE)"
  !insertmacro MUI_INSTALLOPTIONS_WRITE "PythonPage.ini" "Field 2" "State" "$PythonDir"
  !insertmacro MUI_INSTALLOPTIONS_DISPLAY "PythonPage.ini"
FunctionEnd

; Check on startup
Function .onInit
  !insertmacro MUI_INSTALLOPTIONS_EXTRACT "PythonPage.ini"
FunctionEnd

;--------------------------------

; Pages

!insertmacro MUI_PAGE_LICENSE "..\gpl.txt"
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
Page custom PythonPage ;Python custom page
!insertmacro MUI_PAGE_INSTFILES

!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
;Page license
;Page components
;Page custom PythonPage
;Page directory
;Page instfiles

;UninstPage uninstConfirm
;UninstPage instfiles

!insertmacro MUI_LANGUAGE "English"
LangString TEXT_PYTHON_TITLE ${LANG_ENGLISH} "Python Extension Installation"
LangString TEXT_PYTHON_SUBTITLE ${LANG_ENGLISH} "Please locate the directory where Python is installed."

Section "TexGen (required)" ;No components page, name is not important
	
  SectionIn RO

;--------------------------------

; The stuff to install
; Check if 64bit OS
  ${If} ${RunningX64}
	 DetailPrint "Running on x64"
  ${Else}
	 Abort "Not x64 system - please install 32 bit version"
  ${EndIf}

  SetOutPath $INSTDIR

  File TexGenGUI.exe
  File TexGenCore.dll
  File ..\GUI\TexGen.xrc
  File msvcp140.dll
  File vcruntime140.dll
  ;File ..\Docs\TexGen.chm
  File TexGenGUi.exe.manifest
  File Python312\Python312.dll  ; Copy dll anyway in case subversion TexGen is compiled with doesn't match version already installed
  File ..\OctreeRefinement\libsc.dll
  File ..\OctreeRefinement\libp4est.dll
  File ..\OctreeRefinement\zlib1.dll

!insertmacro MUI_INSTALLOPTIONS_READ $PythonDir "PythonPage.ini" "Field 2" "State"
  ;Call CopyPythonDLLIfNeeded

  ;SetOutPath $PythonDir\Lib\site-packages\TexGen
  SetOutPath $INSTDIR\Python\libxtra\TexGen
  
  File ..\Python\__init__.py
  File ..\Python\Lib\Abaqus.py
  File ..\Python\Lib\Ansys.py
  File ..\Python\Lib\Completer.py
  File ..\Python\Lib\FlowTex.py
  File ..\Python\Lib\GridFile.py
  File ..\Python\Lib\WiseTex.py
  File ..\Python\Lib\dataHandling.py
  File ..\Python\Lib\dataHandlingInPlane.py
  File ..\Python\Lib\effectiveMatPropRVE.py
  File ..\Python\Lib\WeavePattern.py
  File ..\Python\Core.py
  File ..\Python\Renderer.py
  File ..\Python\Export.py
  File _Core.pyd
  File _Renderer.pyd
  File _Export.pyd
  File TexGenCore.dll
  File OpenCascade7.7.0\TKBO.DLL
  File OpenCascade7.7.0\TKBOOL.DLL
  File OpenCascade7.7.0\TKBREP.DLL
  File OpenCascade7.7.0\TKERNEL.DLL
  File OpenCascade7.7.0\TKG2D.DLL
  File OpenCascade7.7.0\TKG3D.DLL
  File OpenCascade7.7.0\TKGEOMALGO.DLL
  File OpenCascade7.7.0\TKGEOMBASE.DLL
  File OpenCascade7.7.0\TKIGES.DLL
  File OpenCascade7.7.0\TKMATH.DLL
  File OpenCascade7.7.0\TKOFFSET.DLL
  File OpenCascade7.7.0\TKPRIM.DLL
  File OpenCascade7.7.0\TKSHHEALING.DLL
  File OpenCascade7.7.0\TKSTEP.DLL
  File OpenCascade7.7.0\TKSTEP209.DLL
  File OpenCascade7.7.0\TKSTEPATTR.DLL
  File OpenCascade7.7.0\TKSTEPBASE.DLL
  File OpenCascade7.7.0\TKTOPALGO.DLL
  File OpenCascade7.7.0\TKXSBASE.DLL 
 ; File msvcp71.dll
 ; File msvcr71.dll

  SetOutPath $INSTDIR\Scripts

  File ..\Python\Scripts\2dweave.py
  File ..\Python\Scripts\3dweave.py
  File ..\Python\Scripts\cotton.py
  File ..\Python\Scripts\example.py
  File ..\Python\Scripts\polyester.py
  File ..\Python\Scripts\LayeredTextile.py
  File ..\Python\Scripts\LayeredTextile2.py

  SetOutPath $INSTDIR\Data

  File ..\Data\2dweave.tg3
  File ..\Data\3dweave.tg3
  File ..\Data\cotton.tg3
  File ..\Data\polyester.tg3
  
  SetOutPath $INSTDIR\Utilities
  
  File CFXImportVTK.exe
  File chamis_model_final.for

; Set to 64-bit registry view
  SetRegView 64
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\TexGen "Install_Dir" "$INSTDIR"
  WriteRegStr HKLM SOFTWARE\TexGen "Python_Dir" "$PythonDir"

  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\TexGen" "DisplayName" "TexGen"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\TexGen" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\TexGen" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\TexGen" "NoRepair" 1
  WriteUninstaller "uninstall.exe"
  ; Set back to previous registry view
  SetRegView LastUsed

SectionEnd ; end the section

; Optional section (can be disabled by the user)
Section "Start Menu Shortcuts"
  SetOutPath $INSTDIR   ; Sets the working directory for the links

  CreateDirectory "$SMPROGRAMS\TexGen"
  CreateShortCut "$SMPROGRAMS\TexGen\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  CreateShortCut "$SMPROGRAMS\TexGen\TexGen.lnk" "$INSTDIR\TexGenGUI.exe" "" "$INSTDIR\TexGenGUI.exe" 0
  CreateShortCut "$SMPROGRAMS\TexGen\TexGen API Documentation.lnk" "$INSTDIR\TexGen.chm" "" "$INSTDIR\TexGen.chm" 0

SectionEnd
  
;--------------------------------

; Uninstaller

Section "Uninstall"
  
  ; Set to 64-bit registry view
  SetRegView 64
  ; Find out where the TexGen extensions where installed
  ReadRegStr $PythonDir HKLM SOFTWARE\TexGen "Python_Dir"
  ReadRegStr $INSTDIR HKLM SOFTWARE\TexGen "Install_Dir"

  ; Remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\TexGen"
  DeleteRegKey HKLM "SOFTWARE\TexGen"
  
  ; Set back to previous registry view
  SetRegView LastUsed

  ; Remove files and uninstaller
  Delete $INSTDIR\Python\libxtra\TexGen\__init__.py
  Delete $INSTDIR\Python\libxtra\TexGen\Abaqus.py
  Delete $INSTDIR\Python\libxtra\TexGen\Ansys.py
  Delete $INSTDIR\Python\libxtra\TexGen\Completer.py
  Delete $INSTDIR\Python\libxtra\TexGen\FlowTex.py
  Delete $INSTDIR\Python\libxtra\TexGen\GridFile.py
  Delete $INSTDIR\Python\libxtra\TexGen\WiseTex.py
  Delete $INSTDIR\Python\libxtra\TexGen\Core.py
  Delete $INSTDIR\Python\libxtra\TexGen\Renderer.py
  Delete $INSTDIR\Python\libxtra\TexGen\Export.py
  Delete $INSTDIR\Python\libxtra\TexGen\dataHandling.py
  Delete $INSTDIR\Python\libxtra\TexGen\dataHandlingInPlane.py
  Delete $INSTDIR\Python\libxtra\TexGen\effectiveMatPropRVE.py
  Delete $INSTDIR\Python\libxtra\TexGen\WeavePattern.py

  Delete $INSTDIR\Python\libxtra\TexGen\_Core.pyd
  Delete $INSTDIR\Python\libxtra\TexGen\_Renderer.pyd
  Delete $INSTDIR\Python\libxtra\TexGen\_Export.pyd
  Delete $INSTDIR\Python\libxtra\TexGen\TexGenCore.dll

  Delete $INSTDIR\Python\libxtra\TexGen\TKBO.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKBOOL.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKBREP.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKERNEL.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKG2D.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKG3D.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKGEOMALGO.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKGEOMBASE.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKIGES.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKMATH.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKOFFSET.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKPRIM.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKSHHEALING.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKSTEP.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKSTEP209.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKSTEPATTR.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKSTEPBASE.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKTOPALGO.DLL
  Delete $INSTDIR\Python\libxtra\TexGen\TKXSBASE.DLL
  ;Delete $PythonDir\Lib\site-packages\TexGen\msvcp71.dll
  ;Delete $PythonDir\Lib\site-packages\TexGen\msvcr71.dll

  Delete $INSTDIR\TexGenGUI.exe
  Delete $INSTDIR\TexGenCore.dll
  Delete $INSTDIR\TexGen.xrc
  Delete $INSTDIR\msvcp140.dll
  Delete $INSTDIR\vcruntime140.dll
  Delete $INSTDIR\TexGenGUI.exe.manifest
  ;Delete $INSTDIR\TexGen.chm
  Delete $INSTDIR\Python312.dll
  Delete $INSTDIR\libsc.dll
  Delete $INSTDIR\libp4est.dll
  Delete $INSTDIR\zlib1.dll

  Delete $INSTDIR\Scripts\2dweave.py
  Delete $INSTDIR\Scripts\3dweave.py
  Delete $INSTDIR\Scripts\cotton.py
  Delete $INSTDIR\Scripts\example.py
  Delete $INSTDIR\Scripts\polyester.py
  Delete $INSTDIR\Scripts\LayeredTextile.py
  Delete $INSTDIR\Scripts\LayeredTextile2.py
  
  Delete $INSTDIR\Data\2dweave.tg3
  Delete $INSTDIR\Data\3dweave.tg3
  Delete $INSTDIR\Data\cotton.tg3
  Delete $INSTDIR\Data\polyester.tg3
  
  Delete $INSTDIR\Utilities\CFXImportVTK.exe
  Delete $INSTDIR\Utilities\chamis_model_final.for

  Delete $INSTDIR\uninstall.exe
  Delete $INSTDIR\${PYTHONDLL}

  ; Remove shortcuts, if any
  SetShellVarContext all  ; Doesn't delete shortcuts in Windows 7 without this
  Delete "$SMPROGRAMS\TexGen\*.*"

  ; Remove directories used
  RMDir "$SMPROGRAMS\TexGen"
  ;RMDir "$PythonDir\Lib\site-packages\TexGen"
  RMDir "$INSTDIR\Python\libxtra\TexGen"
  RMDir "$INSTDIR\Python\libxtra"
  RMDir "$INSTDIR\Python"
  RMDir "$INSTDIR\Scripts"
  RMDir "$INSTDIR\Data"
  RMDir "$INSTDIR\Utilities"
  RMDir "$INSTDIR"

SectionEnd
