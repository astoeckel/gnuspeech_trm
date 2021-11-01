# Standalone Gnuspeech Tube Resonance Model (TRM)

This repository contains a standalone version of the Tube Resonance Model (TRM)
used in the Gnuspeech articulatory speech synthesizer. The code comes with a
thin C wrapper that allows to control all aspects of the TRM in
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

In case the Makefile does not work for you, you can compile the executable by
running
```sh
c++ -O3 -s gnuspeech_trm/*.cpp gnuspeech_trm/trm/*.cpp -o build/gnuspeech_trm_main
```
and the library by running
```sh
c++ -O3 -s -fPIC -shared gnuspeech_trm/trm/*.cpp -o build/libgnuspeech_trm.so
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

Please have a look at the code at the end of `gnuspeech_trm/gnuspeech_trm.py`
for the Python wrapper and at `gnuspeech_trm/gnuspeech_trm_main.cpp` for
the for the C wrapper. More details can as well be found in the C header
`gnuspeech_trm/gnuspeech_trm_main.h`.

## License

The code in the `trm` subfolder is licensed under the GPLv3 by various authors.
The C and Python wrappers are (C) 2018 Andreas Stöckel, licensed under the GPLv3.

