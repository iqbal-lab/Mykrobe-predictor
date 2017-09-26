# -*- mode: python -*-

block_cipher = None


a = Analysis(['../mykrobe/mykrobe_predictor.py'],
             pathex=['../../Mykrobe-predictor'],
             binaries=[],
             datas=[],
             hiddenimports=['mykrobe', 'mykatlas'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

a.binaries += Tree('../mccortex/bin/', prefix='.')
a.datas += Tree('../mykrobe/data', prefix='mykrobe/data')

pyz = PYZ(a.pure, a.zipped_data,
          cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='mykrobe_predictor',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='mykrobe_predictor')