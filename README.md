# Cryptanalysis of primitives based on split-and-lookups

## Content

This repository contains implementations of attacks against Tip5 and algorithms to evaluate the linear properties of S-boxes based on split-and-lookups.

* Tip5.ipynb contains a sage implementation of Tip5, an attack against a 4-round reduced version of Tip5 and an integral attack against a full (small field version) of Tip5.
* Split-and-lookup.ipynb contains some tools to evaluate the linear approximations and linear correlation of a split-and-lookup S-box.

## Requirements

You need to have sagemath 9.5 and jupyter notebook installed.

## Running the program

To run the program, run

```
sage -n jupyter
```

in the terminal.

