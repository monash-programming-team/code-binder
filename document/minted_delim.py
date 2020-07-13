# Extract a portion of a file between two given delimiters
#
# Usage:
#
#  python3 minted_delim.py <start> <end>
#
# Source: https://tex.stackexchange.com/questions/130738/formatting-code-fragments-extracted-from-file-with-minted

import sys
codefile = sys.argv[1]
delim1 = sys.argv[2]
delim2 = sys.argv[3]
f = open(codefile)
code = f.read()
f.close()
f = open('snippet.tmp', 'w')  # Save snippet to temp file
try:
    f.write(code.split(delim1, 1)[1].split(delim2, 1)[0])
except:
    f.write('Error trying to split between delimeters "{:s}" and "{:s}"'.format(delim1,delim2))
f.close()
