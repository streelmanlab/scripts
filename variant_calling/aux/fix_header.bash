for file in ~/scratch/ts_ms/combined/*.vcf; do
	echo $file
	grep -v "^##GATKCommandLine" $file > tmp.vcf
	mv tmp.vcf $file
done
#rm ~/scratch/ts_ms/xl_08_27_20/*.idx
