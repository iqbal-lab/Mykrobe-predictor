cd code_generator
python main.py
cd ..
cp code_generator/rendered/include/predictor/core/* include/predictor/core/
cp code_generator/rendered/include/predictor/staph/* include/predictor/staph/
cp code_generator/rendered/src/predictor/core/* src/predictor/core/
cp code_generator/rendered/src/predictor/staph/* src/predictor/staph/