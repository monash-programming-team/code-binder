ASTYLE="./astyle/build/gcc/bin/astyle"

${ASTYLE} --style=google --indent=spaces=2 --indent-namespaces --indent-col1-comments --indent-classes --break-one-line-headers ${1};
