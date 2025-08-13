# Instructions

* ColabFold v1.5.3 (via [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) installation) was run

## MSA Generation

* The exact MMSeqs2 version used was `71dd32ec43e3ac4dabf111bbc4b124f1c66a85f1`

```bash
THREADS="144"
colabfold_search --mmseqs $MMSEQS_PATH/mmseqs --threads $THREADS bphage_and_others_prodigal-gv.faa  $DB_PATH basler_msas
```
* These MSAs were then batched into batches of 500 for folding parallelism across many GPUs i.e.

```bash
cd basler_msas

total_files=$(find . -name "*.a3m" | wc -l)

batch_size=500
total_batches=$(( (total_files + batch_size - 1) / batch_size ))

# Create directories for the batches
for ((i=1; i<=total_batches; i++)); do
    mkdir -p "batch$i"
done

# Move .a3m files to batches
counter=0
batch_number=1

for file in *.a3m; do
    ((counter++))

    # Move to the next batch if the limit is reached
    if [ $counter -gt $batch_size ]; then
        counter=1
	# to tar up the batches
	tar -czf "batch${batch_number}.tar.gz" "batch${batch_number}"  && rm -rf "batch${batch_number}"
        ((batch_number++))
    fi

    # Move the file to the current batch
    mv "$file" "batch$batch_number/"
done
```

## Model Generation

* 3 models were generated for all proteins using 3 recycles (default) and no amber relaxation e.g. 

```
input="batch1"
colabfold_batch  --num-models 3 basler_msas/$input  basler_predictions/$input
```

* Some proteins (the longest ones >3000AA)  didn't fold as the GPUs used (A100/MI250x) didn't have enough VRAM

## Processing

* The script to choose the highest ranking structure, rename and calculate some overall stats (mean pLDDT, protein lengths)

```bash
python rename_process_tophits.py -i bphage_and_others_prodigal-gv.faa -o basler_output_renamed -p basler_predictions/
```


