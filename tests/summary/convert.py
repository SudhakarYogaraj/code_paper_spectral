#!/usr/bin/env python

import os
import sys
from jinja2 import Environment, FileSystemLoader
import yaml

with open(sys.argv[1], 'r') as f:
    data = yaml.load(f)

# Change template
env = Environment(
    block_start_string = '((*',
    block_end_string = '*))',
    variable_start_string = '(((',
    variable_end_string = ')))',
    loader = FileSystemLoader('./templates'))

tex = env.get_template(sys.argv[2]).render(data=data)

with open(sys.argv[2], 'w') as f:
    f.write(tex)
