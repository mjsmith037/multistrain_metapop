# Multi-strain disease dynamics on metapopulation networks

Public code repository associated with the publication [Multi-strain disease dynamics on metapopulation networks](https://www.authorea.com/users/251349/articles/339827-multi-strain-disease-dynamics-on-metapopulation-networks).

This repository is a work in progress; please be forgiving of the current lack of documentation and annotation.

To generate each of the figures from the main text, run the corresponding `r` code in each folder within [/figures](./tree/master/figures). Many of the simulations run by the figure scripts are saved in [/data](./tree/master/data) to reduce re-run time.

Core [code](./tree/master/code) includes:

`multipop_mantis.jl` and `state_space_exploration.jl`, the latter being a wrapper around the former to run over wide parameter ranges.
