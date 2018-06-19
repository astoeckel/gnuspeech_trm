# Standalone Gnuspeech Tube Resonance Model (TRM)

This repository contains a standalone version of the Tube Resonance Model (TRM)
used in the Gnuspeech articulatory speech synthesizer. The code comes with a
thin C wrapper library that allows to control all aspectes of the TRM in
realtime. Furthermore, a Python wrapper is included that allows to integrate the
TRM into any Python application.

## How to build

Building the code was tested on Linux only. You'll need a C++11 compatible
compiler. The code has no dependencies apart from the C++ standard library.

To compile the library, simply run
```sh
git clone https://github.com/astoeckel/gnuspeech_trm
cd gnuspeech_trm
make
```

To test the code, run
```sh
build/gnuspeech_trm_main | aplay -f FLOAT_LE -c 1 -r 44100
```

To install the corresponding Python wrapper, run
```sh
pip3 install --user -e .
```

Test it by running
```sh
gnuspeech_trm/gnuspeech_trm.py | aplay -f FLOAT_LE -c 1 -r 44100
```

## How to use

Please have a look at the bottom of the `gnuspeech_trm/gnuspeech_trm.py` program
for the Python wrapper and at `gnuspeech_trm/gnuspeech_trm_main.cpp` for the 
for the C wrapper. More details can also be found in the C header
`gnuspeech_trm/gnuspeech_trm_main.h`.

## Licence

The code in the `trm` subfolder is licenced under the GPLv3 by various authors.
The C and Python wrappers are (C) 2018 Andreas Stöckel, licenced under the GPLv3.

