### OpenSMOGmod (v2.0)

**OpenSMOGmod** is a patch for [OpenSMOG](https://github.com/smog-server/OpenSMOG) designed to extend the functionality by adding:-
* **Many-particle functions:** Support for custom potentials involving more than 2 atoms.
* **COM Restraints:** Native support for Centre-of-Mass (COM) restraint potentials.

**OpenSMOGmod (version 2.0)** is an **add-on** that uses class inheritance of OpenSMOG SBM class. This allows adding new features without modifying the original source code, thereby ensuring better compatibility with the ongoing development of the base OpenSMOG tool. 
____
#### Installation instrunctions
Since this is an add-on, ensure that OpenSMOG and OpenMM are installed first:

```bash
pip install OpenSMOG
pip install OpenMM
```

Then, clone this repository and install the patch:

```bash
git clone https://github.com/digvijaylp/OpenSMOGmod.git
cd OpenSMOGmod
pip install .
```

Note: The patch has been tested on **OpenSMOG v1.2** (with **OpenMM 8.1**) but can work with older version **OpenSMOG v1.1.1** (with OpenMM 7.7).
____
#### Usage overview
In addition to OpenSMOG input structure (.gro) and topology (.top & .xml) files, OpenSMOGmod requires an additional XML (.xml) file with "OpenSMOGmod" as the root tag. 
```XML
<OpenSMOGmod>
 <manyparticle>
  ...
  ...
 </manyparticle>
 <com_pull>
  ...
  ...
 </com_pull>
</OpenSMOGmod>
```
   Alternatively, the same can be added to the original OpenSMOG XML (.xml) file using XML Processing Instructions.
```XML
<OpenSMOGforces>
 <contacts>
 ...
 </contacts>
<!-- Mod-lines starts -->
 <?OpenSMOGmod>
  <manyparticle>
   ...
  </manyparticle>
  <com_pull>
   ...
  </com_pull>
 ?> <!-- Mod-lines ends -->
</OpenSMOGmod>
```

1. Instead of importing from the base OpenSMOG, import the SBM class from this module. 
```python
from OpenSMOGmod import SBM
```

2. After creating the SBM object, use loadSystemFiles() (from OpenSMOGmod) method instead of OpenSMOG loadSystem(). 
```python 
sbm.loadSystemFiles(Grofile=grofile.gro,Topfile=topfile.top,\
                   Xmlfile=xmlfile,xml,Modfile=modfile.xml)
```
**Note**: Please refer to sample input files and simulation python script in the [examples/](https://github.com/digvijaylp/OpenSMOGmod/examples) directory.
____
#### License
This program is released under the MIT License. A copy of the license is included in the LICENSE file within this repository.

____
#### Usefull likes
[OpenSMOG](https://github.com/smog-server/OpenSMOG): Software for running SBM simulations using OpenMM. <br>
[SuBMIT](https://github.com/sglabncbs/submit): Toolkit for generating CG-SBM input files for MD on OpenSMOG, OpenSMOGmod or GROMACS. <br>
[Legacy OpenSMOGmod](https://github.com/digvijaylp/OpenSMOGmod_legacyfork): Legacy version with modifications made on fork of OpenSMOG v1.1.1.<br>
