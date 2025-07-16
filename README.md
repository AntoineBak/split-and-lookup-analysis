# Cryptanalysis of primitives based on split-and-lookups

## Content

This repository contains implementations of attacks against Tip5 and algorithms to evaluate the linear properties of S-boxes based on split-and-lookups.

* tip_five.sage contains a sage implementation of Tip5, and an attack against a 4-round reduced version of Tip5.
* split_and_lookup.sage contains a class representing split-and-lookup S-boxes and allowing one to compute their linear approximations and correlations.
* compute_vec.mag contains the code to solve the multivariate system generated in the tip_five.sage program.

## Requirements

You need to have sagemath 9.5 and Magma installed.

## Running the program

To run the program, run

```
sage [filename]
```
or
```
magma [filename]
```

in the terminal, depending on the extension of the file.

