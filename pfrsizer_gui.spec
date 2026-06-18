
# -*- mode: python ; coding: utf-8 -*-
block_cipher = None

a = Analysis(
    ['pfrsizer/gui.py'],
    pathex=[],
    binaries=[],
    datas=[('icon.ico', '.')],   # icon bundle injected here (when icon.ico exists in project root)
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
    # Icon for the EXE (taskbar / explorer). Value injected at build time.
    icon=r'C:/Users/Fritz/Desktop/code-projects/pfr-reactor-sizer/icon.ico',
)
