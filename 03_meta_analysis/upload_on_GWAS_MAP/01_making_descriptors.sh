
# Script for automatic generating descriptors json files for GWASes of med Q 
# PGPs done on SOCCS_controls samples
# 

# Path to fodler where descriptors will be stored

descr_folder='/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/json'
mkdir "$descr_folder"

for igp in {1..21}
do
echo $igp

	txt="{"
	echo $txt > "$descr_folder"/igp"$igp".json

	txt="\"trait_abbreviation\":\"IGP"$igp"\","
	echo $txt >> "$descr_folder"/igp"$igp".json

	txt="\"trait_name\":\"synthetic IgG glycan peak "$igp"\","
	echo $txt >> "$descr_folder"/igp"$igp".json

	cat descriptors_end.txt >> "$descr_folder"/igp"$igp".json

done