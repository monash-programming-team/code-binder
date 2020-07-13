# Monash code binder style for Minted
# Attempt to mimick Gedit default colours
#
# You'll need to put a symlink to this file in your minted styles directory
# called monash.py in order to be able to compile the document.

# For example, mine is symlinked from
#
#   ~/.local/lib/python3.5/site-packages/pygments/styles/monash.py

from pygments.style import Style
from pygments.token import Keyword, Name, Comment, String, Error, \
     Number, Operator, Generic, Literal

class MonashStyle(Style):
    default_style = ""
    styles = {
        Comment:             '#00F',           # Code comments
        Comment.Preproc:     '#a220f0',        # Preprocessor macros
        Comment.PreprocFile: '#ff0cff',        # Include file
        Keyword:             'bold #a52a2a',   # General keywords
        Keyword.Type:        'bold #39794c',   # Types
        Name:                '#000',           # Variable and function names
        Literal:             '#ff0cff',        # String and numeric literals
    }

