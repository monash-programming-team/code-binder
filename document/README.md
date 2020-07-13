# TEAM REFERENCE DOCUMENT #

This document consists of a compilation of (mostly) well-tested code that can be used as a reference document in official ICPC contests.

## Dependencies for building the document ##

The latest document should be stored in the repository, so if you only want to read the document, you don't need to build anything. If you want to modify the document, then you'll need to be able to build it. To build the document, you will need:

* A LaTeX compiler
 * The [minted](https://ctan.org/pkg/minted?lang=en) package
* Python3
 * The [pygments](http://pygments.org/) library (version 2.2.0 or newer)


### Installing Minted ###

Minted is included in standard [TeXLive](https://ctan.org/pkg/texlive) and [MiKTeX](https://ctan.org/pkg/miktex) distributions.

### Installing Pygments ###

The easiest way to install pygments is using [pip](https://packaging.python.org/tutorials/installing-packages/), running

```pip install --user pygments```

If you have an older version of Pygments (type `pygmentize -V` to check,) you should update it with

```pip install --user --upgrade pygments```

### Installing our Custom Style ###

The Monash team reference document uses a custom style for Pygments. You can install it by running the included setup script `setup.sh`, which will create a symbolic link to the style in your Pygments installation.

## Building ##

To build the document, you will need to run LaTeX with the flag `--shell-escape` so that Pygments can be run. For example, you could compile the document with

```pdflatex --shell-escape document.tex```

or with your favourite LaTeX IDE, as long as you have enabled the `--shell-escape` flag.
