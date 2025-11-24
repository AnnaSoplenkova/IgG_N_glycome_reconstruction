# запускать из /mnt/polyomica/projects/GPD

  ./GlycomeQC \
  --project TrainIgG_merge45_2025_QC \
  --data /home/common/projects/IgG_prediction_from_Plasma/data/TwinsLongitudinalIgG_raw_data_corrected.txt \
  --out /mnt/polyomica/projects/GPD/training_set/QC_results \
  --norm mq \
  --tmp /mnt/polyomica/projects/GPD/training_set/QC_results/QC_output_tmp \
  --harm GP4-GP5
