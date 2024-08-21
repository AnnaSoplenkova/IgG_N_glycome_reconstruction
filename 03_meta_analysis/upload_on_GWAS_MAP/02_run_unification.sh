descr_folder='/home/soplenkova/polyomica/projects/GPD/result_ma/IgG_OLS/json'
uni_folder='/home/soplenkova/polyomica/projects/GPD/result_ma/IgG_OLS/unification_out'


for gp in {1..21}
do
echo $gp

mkdir "$uni_folder"/igp"$gp"

python3 ${PYTHONPATH}/gwas_v2/unifier/run_uni_qc_rep.py \
--gwas-path=/home/soplenkova/polyomica/projects/GPD/result_ma/IgG_OLS/10k_IGP_"$gp".csv \
--mapping-path=mapping.json \
--descriptors-path="$descr_folder"/igp"$gp".json \
--output-dir="$uni_folder"/igp"$gp"/ \
--output-file=igp"$gp" \
--qc-report

done