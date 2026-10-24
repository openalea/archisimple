# ArchiSimple

ArchiSimple is a parsimonious model of root system architecture (Pagès et al., 2014): each apex elongates in
proportion to its diameter, for a growth duration set by its diameter (determinate laterals); laterals form at a regular
distance behind the apices, their diameters drawn from their bearer's; older roots thicken by the pipe model; and the
root system's growth is the potential one times a satisfaction coefficient, the structural mass offered over the mass
the potential growth demands.

This package implements it as a [metafspm](https://github.com/openalea/metafspm) StructuralComponent: the root system
lives in metafspm's multiscale plant graph (MPG), its steps vectorised on the DataStructure's arrays, so that plants of
a population grow together in one graph. It is the architecture that [RhizoDep](https://github.com/openalea/rhizodep)'s
root growth specialises (growth limited by resources, root hairs).

Pagès L., Bécel C., Boukcim H., Moreau D., Nguyen C., Voisin A.-S. (2014). Calibration and evaluation of ArchiSimple,
a simple model of root system architecture. *Ecological Modelling* 290, 76–84.

## Installation

metafspm is not installable from `openalea3` in this version yet: clone the repositories side by side and install them
from source.

```bash
mkdir fspm && cd fspm
git clone -b data_structure_api https://github.com/openalea/metafspm.git metafspm
git clone https://github.com/openalea/archisimple.git archisimple

mamba env create -f archisimple/conda/environment.yml
mamba activate archisimple
pip install -e metafspm
pip install -e archisimple
```
