"""
Build a standalone Windows .exe for the PFR Reactor Sizer GUI.

Usage (on Windows with PyInstaller installed):
    python build_exe.py

This produces a single-file executable (dist/pfrsizer-gui.exe or similar)
that can be distributed to other Windows PCs (no Python installation needed).

Requirements for building:
    pip install pyinstaller

Notes:
- The build includes numpy, scipy, matplotlib, pubchempy, etc.
- Matplotlib and scipy often require several --hidden-import flags.
- The resulting EXE will be fairly large (~150-250 MB) because of the scientific stack.
"""

import os
import sys
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
DIST_DIR = PROJECT_ROOT / "dist"

SPEC = """
# -*- mode: python ; coding: utf-8 -*-
block_cipher = None

a = Analysis(
    ['pfrsizer/gui.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=[
        'numpy',
        'scipy',
        'scipy.integrate',
        'scipy._lib',
        'matplotlib',
        'matplotlib.backends.backend_tkagg',
        'pubchempy',
        'requests',
        'tkinter',
        'tkinter.ttk',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='PFR_Reactor_Sizer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,                 # Windowed app (no console)
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,                     # Add 'icon.ico' here if you have one
)
"""

def main():
    print("=== Building PFR Reactor Sizer Windows EXE ===")
    print("Project root:", PROJECT_ROOT)

    # Ensure pyinstaller is available
    try:
        import PyInstaller
    except ImportError:
        print("PyInstaller not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])

    # Write a .spec file for full control
    spec_path = PROJECT_ROOT / "pfrsizer_gui.spec"
    spec_path.write_text(SPEC, encoding="utf-8")
    print("Wrote spec file:", spec_path)

    # Clean previous builds
    for p in [PROJECT_ROOT / "build", DIST_DIR]:
        if p.exists():
            import shutil
            shutil.rmtree(p, ignore_errors=True)

    cmd = [
        sys.executable, "-m", "PyInstaller",
        "--clean",
        "--noconfirm",
        str(spec_path),
    ]
    print("Running:", " ".join(cmd))
    subprocess.check_call(cmd, cwd=PROJECT_ROOT)

    print("\n=== Build complete ===")
    exe_path = DIST_DIR / "PFR_Reactor_Sizer.exe"
    if exe_path.exists():
        print("Executable:", exe_path)
        print("You can copy this .exe to any Windows PC and run it directly.")
    else:
        print("Check the dist/ folder for the output executable.")


if __name__ == "__main__":
    main()
