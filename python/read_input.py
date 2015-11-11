import os
import sys

# Keywords to substitute
keywords = ["solution", "second order drift", "potential", "dimensions", "diffusion"]

# Only 1 input file
input_file = sys.argv[1]
output_file = sys.argv[2]

# Copy build_problem_init.py
os.system("cp build_problem_init.py " + output_file)

for keyword in keywords:

    # Open files and read lines
    input_lines = open(input_file).readlines()
    output_lines = open(output_file).readlines()

    # Open output file in 'write' mode
    foutput = open(output_file, 'w')

    # Boolean to know if we are within one piece
    in_input_text = False

    # Input corresponding to keyword
    user_input = ""

    # Read in input file
    for line in input_lines:
        if line.startswith("# user input : " + keyword):
            in_input_text = True
        if in_input_text:
            user_input += line
        if line.startswith("# end"):
            in_input_text = False

    # Write in output file
    for line in output_lines:
        if line.startswith("# user input : " + keyword):
            in_input_text = True
        if not in_input_text:
            foutput.write(line)
        if line.startswith("# end") & in_input_text:
            in_input_text = False
            foutput.write(user_input)

    # Close output file
    foutput.close()
