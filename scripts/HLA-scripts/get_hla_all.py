#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7
import os
import csv

input_root = 'processed'
output_file = 'aggregated_hla_genotypes.txt'

hla_results = []

for root, dirs, files in os.walk(input_root):
    for file in files:
        if file.endswith('_result.tsv'):
            tsv_path = os.path.join(root, file)
            sample_name = os.path.basename(os.path.dirname(root)) + '_1.bed'

            with open(tsv_path) as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    hla_fields = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
                    hla_values = ['HLA-' + row[field] for field in hla_fields if row[field]]
                    hla_string = ','.join(hla_values)
                    hla_results.append((sample_name, hla_string))
                    break  # Only process the first row

# Write results to file
with open(output_file, 'w') as out:
    out.write('sample\thla\n')
    for sample, hla in hla_results:
        out.write(f'{sample}\t{hla}\n')
