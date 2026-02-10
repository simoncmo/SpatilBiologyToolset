# 2025-05-26 Simon
echo "Backup This Folder To Box"
# Setup
katmai_dir="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis"
box_backup_dir="0_Katmai_backup"
box_name="Box_Dinglab"
n_files_at_once=10

# Run
rclone mkdir "${box_name}":"${box_backup_dir}/$(basename ${katmai_dir})"
rclone copy -P --copy-links --transfers ${n_files_at_once} --exclude="*.fastq*" "${katmai_dir}" "${box_name}":"${box_backup_dir}/$(basename ${katmai_dir})"


# Notes:
# -L/--copy-links: follows symlinks
# -transfers 10: transfer 10 files at once
# --exclude="*.fastq*": Ignore fastq files
