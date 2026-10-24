# ArchiSimple

ArchiSimple (Pagès et al., 2014) is a parsimonious model of root system architecture: the apices' elongation from
their diameter, determinate growth durations, the branching density and the laterals' diameters from their bearer's,
radial growth by the pipe model, and the growth limited by a satisfaction coefficient. This package implements it as a
`metafspm` StructuralComponent, vectorised on the plant's multiscale graph; RhizoDep's root growth specialises it.

```{include} ../README.md
:start-after: "## Installation"
```

## Tutorial

`example/tutorials/tutorial_archisimple.ipynb`: a plant model and its scenario, the MPG's grouping, the satisfaction
coefficient, a population in one MPG, the calibrated species.

## Examples

- `example/species.py`: the calibrated species, by root order and by radius;
- `example/profile_population.py`: a population in one MPG against its plants grown one by one.
