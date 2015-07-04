# Select lines in user input
inputs = ["", "", "", ""]
finput = open('user_input.py', 'r')
lines = finput.readlines()
finput.close()

index_input = 0
in_input = False

for line in lines:
    if line.startswith("# USER INPUT"):
        in_input = True
    if in_input == 1:
        inputs[index_input] += line
    if line.startswith("# END"):
        in_input = False
        index_input += 1

# Read lines of build_problem.py
fbuild = open('build_problem.py', 'r')
lines = fbuild.readlines()
fbuild.close()

# Write input in build_problem.py
fbuild = open('build_problem.py', 'w')

index_input = 0
in_input = False

for line in lines:
    if line.startswith('# USER INPUT'):
        in_input = True
    if not in_input:
        fbuild.write(line)
    if line.startswith('# END OF USER INPUT'):
        fbuild.write(inputs[index_input])
        in_input = False
        index_input += 1
