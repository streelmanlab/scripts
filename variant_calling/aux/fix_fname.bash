for f in ~/scratch/ts_ms/xl_08_28_20/*.vcf;do
	contig=$(tail -n 1 $f | cut -f 1)
	mv $f ~/scratch/ts_ms/xl_08_28_20/xl_${contig}.vcf
	mv "$f.idx" ~/scratch/ts_ms/xl_08_28_20/xl_${contig}.vcf.idx
	#echo $f
	#bn="$(basename -- $FILE)"
done
