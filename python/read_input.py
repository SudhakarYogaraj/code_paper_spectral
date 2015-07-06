def set_in_file(keyword):

    # Open files and read lines
    input_lines = open('user_input.py').readlines()
    output_lines = open('build_problem.py').readlines()

    # Open output file in 'write' mode
    foutput = open('build_problem.py', 'w')

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

keywords = ["solution", "second order drift", "potential", "dimensions"]

for keyword in keywords:
    set_in_file(keyword)
