#!/bin/bash

# =====================================
# To download multiple files listed as curl commands in a file
# =====================================

# File that contains a list of curl commands, one per line
URL_FILE="wgs_top_50_urls.sh"

# Directory where downloaded files will be stored
TARGET_DIR="raw_data_wgs"

# ----------------------------------------------------------
# STEP 1: Check if the input file exists
# -----------------------------------------------------------
if [[ ! -f "$URL_FILE" ]]; then
  echo "ERROR: File '$URL_FILE' not found!"
  exit 1                                           
 # Exit the script if the file is not found
  else
  echo "File '$URL_FILE' found. Starting downloads"
fi

# -----------------------------------------------------------
# STEP 2: Create the target directory if it doesn't exist already
# ------------------------------------------------------------
mkdir -p "$TARGET_DIR"  
# -p ensures no error if the directory already exists

# ---------------------------------------------------------------
# STEP 3: Loop through each line of the file
# ---------------------------------------------------------------
# while IFS= Read the file line-by-line, keeping the line exactly as it is without breaking it at spaces    
# read -r says: donot interpret backslashes as escape characters
while IFS= read -r line; do
 # Skip empty lines or lines starting with #
  [[ -z "$line" || "$line" =~ ^# ]] && continue

  # -----------------------------------------------------------------
  # STEP 4: Extract the filename from the curl command using grep
  # look behind regex (?<=-o), looks for the string after '-o ' and [^ ]+ match all non-space      characters after -o
  # -----------------------------------------------------------------
  OUTPUT_FILE=$(echo "$line" | grep -oP '(?<=-o )[^ ]+')

  # --------------------------------------------------------------------
  # STEP 5: Set the full new path for the output file
  # -------------------------------------------------------------------
 BASENAME=$(basename "$OUTPUT_FILE")  # skips the path and keeps the file name
 NEW_OUTPUT="$TARGET_DIR/${BASENAME}" # sets the path

  # -------------------------------------------------------------------
  # STEP 6: the original output file in the curl command is replaced with the new one which has the path.
  # This ensures all files are saved in the target directory
  # ----------------------------------------------------------------------
  MODIFIED_CURL=$(echo "$line" | sed "s|-o $OUTPUT_FILE|-o $NEW_OUTPUT|")

  # ------------------------------------------------------------------------
  # STEP 7: Print the destination file and run the curl command
  # -------------------------------------------------------------------------
  echo "Downloading to: $NEW_OUTPUT"
  eval "$MODIFIED_CURL"      # Run the modified command using eval

done < "$URL_FILE"       # Feed the file line by line to the while loop

# ------------------------------------------------------------------------
# STEP 8: Done message
# ------------------------------------------------------------------------
echo "All files downloaded to '$TARGET_DIR/'"


