# geltools

## Compile

```
mkdir build
cd build
cmake ..
make
```

The geltools executable will appear in build/src/geltools. To run the test, one can run `make test`.

## Example

```
cd build

# when the imputed data is not provided with dosage. (use --no_dosage (anything))
./src/geltools --input ../tests/winni_regression/test.vcf.gz --truth ../tests/winni_regression/truth.vcf.gz --freq ../tests/winni_regression/test_geltools.freq  --mode r2 --no_dosage x 

# when the imputed data is provided with dosage. e.g. produced by IMPUTE5 or Beagle5

./src/geltools --mode r2 --input ../tests/beagle_out_samll.vcf --truth ../tests/beagle_out2.vcf --freq ../tests/beagle_out2.af

```
