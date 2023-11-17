This forked version is for PUEO, superceding the previous pueo branch


# Models of Antarctica in ROOT
For analysis and simulation of the ANITA experiment.
Although, much of the utilities here could apply to any Antarctic based radio neutrino experiment.

## The backstory

Analysis of ANITA data and simulation of ANITA both require models of Antarctica.
The existence of this library is an effort to deduplicate code during a comprehensive refactor of icemc (the simulation code).
Since the analysis code was in a better state, this repository begins is a fork of the anitaEventCorrelator library, with the ANITA analysis dependent things removed.
The history of those classes in this repository was `git-filter-branch`-ed away, so beyond a certain point the commit history here might not make much sense.
However, if you really need to jump back in time the history of these classes at the time of the fork is still available in the anitaEventCorrelator library.

## The future
Models from icemc not already  present in this library will also make their way here.
