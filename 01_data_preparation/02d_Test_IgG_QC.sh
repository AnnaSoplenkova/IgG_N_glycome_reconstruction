# запускать из /mnt/polyomica/projects/GPD
./GlycomeQC \
  --project TestIgG_merge45_2025_QC \
  --data /mnt/polyomica/projects/GPD/test_set/SYSCID_CEDAR_IgG_for_QC.csv \
  --out /mnt/polyomica/projects/GPD/test_set/QC_results \
  --norm mq \
  --tmp /mnt/polyomica/projects/GPD/test_set/QC_results/QC_output_tmp \
  --harm GP2-GP3,GP5-GP6,GP13-GP14
