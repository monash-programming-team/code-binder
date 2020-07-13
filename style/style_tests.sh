for file in $(find ../official -name \*.cpp); do
  ./style_cpp.sh ${file};
  if [ -f ${file}.orig ]; then
    diff ${file} ${file}.orig || exit 1
  fi
done
