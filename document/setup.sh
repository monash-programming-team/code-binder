# Set up the custom pygments style required for building the document
# 
# Usage: setup.sh

# ----------------------------------------------------------------------
# 										Locate the pygments installation
# ----------------------------------------------------------------------

pygments_file=`python3 -c "import pygments; print(pygments.__file__)"`

if ! grep -i "__init__.py" <(echo $pygments_file) >/dev/null; then
	echo "Could not load required Python3 module <pygments>"
	echo "Try installing it via 'pip install --user pygments'"
	exit
fi

pygments_loc=$(dirname $pygments_file)
echo -e "\e[32;1mGood\e[0m: Found <pygments> module at ${pygments_loc}"

# ----------------------------------------------------------------------
#												"Install" our custom style
# ----------------------------------------------------------------------

## This line may require sudo depending on your system setup
ln -s -f $(pwd)/style.py ${pygments_loc}/styles/monash.py || exit
echo -e "\e[32;1mGood\e[0m: Created symlink to custom pygments style"

