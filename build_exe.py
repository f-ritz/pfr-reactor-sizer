"""
Build a standalone Windows .exe for the PFR Reactor Sizer GUI.

Usage (on Windows with PyInstaller installed):
    python build_exe.py

This produces a single-file executable (dist/pfrsizer-gui.exe or similar)
that can be distributed to other Windows PCs (no Python installation needed).

Requirements for building:
    pip install pyinstaller

Notes:
- The build includes numpy, scipy, matplotlib, pubchempy, requests, etc.
- Matplotlib and scipy often require several --hidden-import flags.
- The resulting EXE will be fairly large (~150-250 MB) because of the scientific stack.
- To set a custom icon: drop icon.ico (256x256 or multiple sizes recommended) in the project root before running build_exe.py.
  It is used for the EXE icon (file properties) and also bundled so the Tk GUI can set the runtime window + taskbar icon via iconbitmap.
  If the taskbar icon doesn't update after rebuild, restart explorer.exe or clear the Windows icon cache.
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
    datas={datas},   # icon bundle injected here (when icon.ico exists in project root)
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
    hooksconfig={{}},
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
    # Icon for the EXE (taskbar / explorer). Value injected at build time.
    icon={icon},
)
"""

def main():
    print("=== Building PFR Reactor Sizer Windows EXE ===")
    print("Project root:", PROJECT_ROOT)

    # Check if the old EXE is still running
    exe_path = DIST_DIR / "PFR_Reactor_Sizer.exe"
    if exe_path.exists():
        try:
            # Try to rename it — if it fails, it's likely still running
            temp_path = exe_path.with_suffix(".exe.old")
            exe_path.rename(temp_path)
            temp_path.unlink()
        except PermissionError:
            print("WARNING: The previous EXE appears to still be running.")
            print("Please close PFR_Reactor_Sizer.exe before rebuilding.")
            print("Continuing anyway (will retry deletions)...")

    # Ensure pyinstaller is available
    try:
        import PyInstaller
    except ImportError:
        print("PyInstaller not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])

    # Compute icon path at build time (in normal Python, where __file__ works)
    # We embed a concrete value so the generated .spec never relies on __file__ or PROJECT_ROOT.
    icon_path = PROJECT_ROOT / "icon.ico"
    if icon_path.exists():
        # Use raw string with forward slashes (works reliably in the spec)
        icon_literal = f"r'{icon_path.as_posix()}'"
        datas_literal = "[('icon.ico', '.')]"
    else:
        icon_literal = "None"
        datas_literal = "[]"
        print("No icon.ico found in project root — building without custom icon.")

    # Write a .spec file for full control
    spec_path = PROJECT_ROOT / "pfrsizer_gui.spec"
    if spec_path.exists():
        try:
            spec_path.unlink()
        except Exception:
            pass
    spec_content = SPEC.format(icon=icon_literal, datas=datas_literal)
    spec_path.write_text(spec_content, encoding="utf-8")
    print("Wrote spec file:", spec_path)
    print(f"  Icon setting: {icon_literal}")
    print(f"  Datas: {datas_literal}")

    # Clean previous builds — with retries for locked files
    import shutil
    import time

    def _rmtree_retry(path, max_retries=5, delay=1.0):
        """Remove a directory tree, retrying if files are temporarily locked."""
        for attempt in range(max_retries):
            try:
                if path.exists():
                    shutil.rmtree(path)
                return
            except PermissionError:
                if attempt < max_retries - 1:
                    print(f"  File locked, retrying in {delay}s... (attempt {attempt + 1}/{max_retries})")
                    time.sleep(delay)
                else:
                    # Last resort: try to delete individual files
                    print(f"  Could not remove {path} entirely, trying individual files...")
                    for item in path.rglob("*"):
                        try:
                            if item.is_file():
                                item.unlink()
                        except PermissionError:
                            print(f"  Warning: could not delete {item} (still in use)")
                    try:
                        shutil.rmtree(path)
                    except PermissionError:
                        print(f"  Warning: some files in {path} could not be removed (still in use)")

    for p in [PROJECT_ROOT / "build", DIST_DIR]:
        _rmtree_retry(p)

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
