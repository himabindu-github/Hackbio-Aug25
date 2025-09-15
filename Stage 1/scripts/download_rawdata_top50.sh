#!/bin/bash

# Input file containing curl commands
URL_FILE="wgs_top_50_urls.sh"
TARGET_DIR="raw_data_wgs"

# Check if file exists
if [[ ! -f "$URL_FILE" ]]; then
  echo " ERROR: File '$COMMANDS_FILE' not found!"
  exit 1
else
  echo "File '$URL_FILE' found. Starting downloads..."
fi

# Create target folder if not exists
mkdir -p "$TARGET_DIR"

# Read and execute each curl command
while IFS= read -r line; do
  # Skip empty lines or comments
  [[ -z "$line" || "$line" =~ ^# ]] && continue

  # Extract the filename after -o
  OUTPUT_FILE=$(echo "$line" | grep -oP '(?<=-o )[^ ]+')

  # Set new output path
  NEW_OUTPUT="$TARGET_DIR/${OUTPUT_FILE##*/}"

  # Replace original output path in command
  MODIFIED_CMD=$(echo "$line" | sed "s|-o $OUTPUT_FILE|-o $NEW_OUTPUT|")

  # Print and execute
  echo "Downloading to: $NEW_OUTPUT"
  eval "$MODIFIED_CMD"

done < "$URL_FILE"

echo "All files downloaded to $TARGET_DIR/"
