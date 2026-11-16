# nwergm-tripartite (development version of nwergm)

This folder is a **development fork of the Stata package `nwergm`** that adds the model for exclusive multi-affiliation (agency–player–club) networks used in the paper: the typed-swap sampler, the cross-layer terms (`tpcoloc`, `tpstarcov`, ...), and degree-preserving constraints. It is **not** the released `nwergm`. The command names are unchanged (`nwergm`), so **do not install it alongside the released package**: use it only through `adopath ++` from this repository, as the replication scripts do (`do-files` add `software/nwergm-tripartite` and `software/nwergm-tripartite/lib/plugins/macos` to the ado path).

**Requires the core package `nwcommands` (3.6.0 or later)**: in Stata, `net install nwcommands, from(https://nwcommands.org/net/) replace`.

The native MCMC plugins are in `lib/plugins/` (macOS, Linux, Windows); the replication runs used the macOS plugin. License: GPL-3 (`LICENSE`).
