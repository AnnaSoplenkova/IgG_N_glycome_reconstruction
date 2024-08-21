
descr_folder='/home/soplenkova/polyomica/projects/GPD/result_ma/IgG_OLS/json'
uni_folder='/home/soplenkova/polyomica/projects/GPD/result_ma/IgG_OLS/unification_out'

for gp in {1..21}
do
echo $gp

python3 ${PYTHONPATH}/gwas_v2/unifier/run_upload.py \
--gwas-path="$uni_folder"/igp"$gp"/igp"$gp"_done.csv

done