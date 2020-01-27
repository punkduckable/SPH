#!/bin/bash

# Initialize totals to 0
total=0
source_subtotal=0
test_subtotal=0


################################################################################
# Find source subtotal

# First, put all .cc and .h files from the source directory in a list.
SOURCE_FILE_LIST=$(find ./src | grep  -E "(.h$|.cc$)")

# Now loop through the files, adding the length of each one to the subtotal
for File in $SOURCE_FILE_LIST
do
  # Get length of file by counting the number of new lines (^)
  line_count=$(grep -c ^ $File)

  # Add number of lines to subtotal
  let source_subtotal=source_subtotal+line_count
done


################################################################################
# Find test subtotal

# First, put all .cc and .h files from the test directory in a list
TEST_FILE_LIST=$(find ./test | grep -E "(.h$|.cc$)")

# Now loop through the files, adding the length of each one to  the subtotal
for File in $TEST_FILE_LIST
do
  # Get length of file by counting the number of new lines (^)
  line_count=$(grep -c ^ $File)

  # Add number of lines to subtotal
  let test_subtotal=test_subtotal+line_count
done


################################################################################
# Report totals

Makefile_subtotal=$(grep -c ^ ./Makefile)
let total=source_subtotal+test_subtotal+Makefile_subtotal

printf "lines in Source:             %s\n" $source_subtotal
printf "lines in Test:               %d\n" $test_subtotal
printf "Total (including Makefile):  %d\n" $total
