# Instructions

# CF v1.5.5

# MSA 

```

colabfold_search --mmseqs $MMSEQS_PATH/mmseqs --threads $THREADS bphage_and_others_prodigal-gv.faa  $DB_PATH basler_msas


```

# Inference

```
colabfold_batch  --num-models 3 basler_msas/$input  basler_predictions/$input

```

* I put all files into `basler_predictions`
* Very few files didn't fold (too long)

# Processing

python rename_process_tophits.py -i bphage_and_others_prodigal-gv.faa -o basler_output_renamed -p basler_predictions/


# Tarballs

* Massive tarball was 31GB, too large to share. Therefore, used `split`

```
split -b 11G basler_output_renamed.tar.gz "basler_output_renamed.tar.gz.part"
```

* To put back together again

```
cat basler_output_renamed.tar.gz.part* > basler_output_renamed.tar.gz
tar -xzf basler_output_renamed.tar.gz
```


