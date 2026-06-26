#!/bin/bash

# Change param in TmatParam files

#initialname="improved"
#newname="David2026PhD"

initialname="elev90"
newname="vertical"


# Loop through all files in current directory
for file in *; do
    # Check if file name contains the initial name
    if [[ "$file" == *"$initialname"* ]]; then
        # Replace the initial name with new name
        new_filename="${file//${initialname}/${newname}}"
        # Rename the file
        mv "$file" "$new_filename"
        echo "Renamed: $file -> $new_filename"
    fi
done