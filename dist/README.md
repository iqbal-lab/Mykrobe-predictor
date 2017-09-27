# Rough procedure for generating the win64 installer (`mykrobe_predictor_installer_win64.exe`)

Working directory: `Mykrobe-predictor`

1) modify `./mccortex/libs/xxHash/xxhsum.c` according to https://github.com/Cyan4973/xxHash/issues/100

2) with cygwin64: `(cd ./mccortex && make)`

3) with cygwin64: `pip install git+https://github.com/Phelimb/atlas`

4) with cygwin64: `pip install pyinstaller`

5) with cygwin64: `(cd dist && pyinstaller --workpath='./pyinstaller_build/binary_cache' --distpath='./pyinstaller_build' mykrobe_predictor_pyinstaller.spec)`

6) install InnoSetup (http://www.jrsoftware.org/isinfo.php)

7) compile the install binary using the InnoSetup GUI and `mykrobe_predictor_win64.iss`

8) win64 installer binary: `./dist/inno_setup_build/mykrobe_predictor_installer_win64.exe`