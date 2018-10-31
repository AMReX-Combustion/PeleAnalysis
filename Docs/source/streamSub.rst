.. highlight:: bash


streamSub
*********

Extract a subset of streamlines from a streamline file, based no a number of criteria (random sampling, physical
location of the streamline seed point, etc). This process will fail to respect the connectivity of the streamlines
so that the result will no represent the bounds of subvolumes in the domain.  Typically, this tool is use to
extract a management number of lines for quick plotting in order to get a feel for the data contained in the
full set (or to conditionally sample a larger set to gather statistics of data on the streamlines).


```
Usage:
```

Example:

