./measureAggregateRsquare --validation truth.gen.gz --imputed test.gen.gz --sample truth.samples --freq test.freq --bin olivierBins.txt --output test
./measureAggregateRsquare --validation truth_one_row.gen.gz --imputed test_one_row.gen.gz --sample truth_one_row.samples --freq test_one_row.freq --bin olivierBins.txt --output test_one_row

../../build/src/geltools --input test.vcf.gz --truth truth.vcf.gz --freq test_geltools.freq  --mode r2 --no_dosage x
cat test_geltools.freq.r2

../../build/src/geltools --input test_one_row.vcf --truth truth_one_row.vcf --freq test_one_row_geltools.freq  --mode r2 --no_dosage x
cat test_one_row_geltools.freq.r2

