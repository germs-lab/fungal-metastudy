# fungal-metastudy

### Prepare file for co-occurence analysis
R script is used to prepare file in a correct format

### Runing co-occurence analysis
Compile
```
c++ get-co-occurrence-table.cpp -o co -fopenmp
```
Run
```
./co -f tomas_table_for_co.unix.csv -t 2 -s 3 -p 0.05 > tomas_result_p0.05.tsv
```

### Filter rho value over 0.6
```
python filter_co_result.py tomas_result_p0.05.tsv 0.6 > tomas_result_p0.05_r0.6.tsv
```

### Group
R script is used to group

