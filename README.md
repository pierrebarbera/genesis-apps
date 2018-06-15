# Applications for the Genesis Library

A collection of apps for use with the [genesis library for phylogenetic placement post-analysis](https://github.com/lczech/genesis).

## Usage:
Simply navigate to the genesis/apps folder and clone this repository into it. Genesis will recognize it as an application subdirectory and build the apps with the next call to `make` (or `make update` if the build directory exists).

## Rules for Contributing
The style of these apps is pretty loose, but some rules should be followed:

- Filenames should be descriptive, dashes should be used as word separators
- apps should have one main purpose (UNIX-style)
- ideally apps should print to stdout such that they can be piped
