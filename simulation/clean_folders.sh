#!/bin/bash
# Define the list of directories
directories=(result collected_result output output_diag)

# Loop over each directory
for dir in "${directories[@]}"; do
  if [ -d "$dir" ]; then
    # Remove only files (not subdirectories) in the directory (non-recursive)
    find "$dir" -maxdepth 1 -type f -exec rm -f {} \;
    echo "Removed files from $dir"
  else
    echo "Directory $dir does not exist"
  fi
done
