# Set the directory containing the files
directory <- "/Users/bopeng/Desktop/Lab/Givaudan/GivaudanRun1_16S"

# List all files in the directory
files <- list.files(directory)

# Define the pattern to match "_S***_" in file names
pattern <- "_S.*?_"

# Loop through each file
for (file in files) {
  # Check if the file name matches the pattern
  if (grepl(pattern, file)) {
    # Extract the part between "S" and "_"
    replacement <- gsub(pattern, "_", file)
    # Rename the file
    file.rename(file.path(directory, file), file.path(directory, replacement))
  }
}



