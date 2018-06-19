#!/usr/bin/env python3

#   gnuspeech_trm -- Standalone gnuspeech articulatory synthesiser
#   Copyright (C) 2018  Andreas Stöckel
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
Python wrapper for gnuspeech_trm. Uses the foreign funtion interface to provide
a class-based wrapper around gnuspeech_trm
'''

import os
import ctypes
import numpy as np


class TRM:
    # Name of the library the TRM class is supposed to load
    _libname = 'libgnuspeech_trm.so'

    # Handle at the native library. This is lazily initialised
    _lib = None

    # Number of sections in the tube model
    _total_sections = 8

    # Number of nose coefficients
    _total_nose_coefs = 5

    # Number of parameters (not configuration options)
    _total_params = 16

    # Type used instead of ctypes._c_void_p, which does not work for some reason
    # (clamps pointers to 32 bot)
    _c_void_p = ctypes.POINTER(ctypes.c_byte)

    # Pointer at a 32-bit floating point number
    _c_float_p = ctypes.POINTER(ctypes.c_float)

    # Pointer at a 64-bit floating point number
    _c_double_p = ctypes.POINTER(ctypes.c_double)

    # Map from configuration option names to the corresponding key
    _conf_map = {
        'output_rate': 10,
        'control_rate': 20,
        'filter_period': 30,
        'volume': 40,
        'waveform': 50,
        'temperature': 60,
        'loss_factor': 70,
        'mouth_coeff': 80,
        'nose_coeff': 90,
        'throat_cutoff': 100,
        'throat_vol': 110,
        'modulation': 120,
        'mix_offset': 130,
        'glottal_pulse_tp': 140,
        'glottal_pulse_tn_min': 150,
        'glottal_pulse_tn_max': 160,
        'breathiness': 170,
        'vocal_tract_length': 180,
        'reference_glottal_pitch': 190,
        'aperture_radius': 200,
        'nose_radius': (1000, _total_nose_coefs),
        'nose_radius_1': 1000,
        'nose_radius_2': 1001,
        'nose_radius_3': 1002,
        'nose_radius_4': 1003,
        'nose_radius_5': 1004,
        'radius_coef': (2000, _total_sections),
        'radius_coef_1': 2000,
        'radius_coef_2': 2001,
        'radius_coef_3': 2002,
        'radius_coef_4': 2003,
        'radius_coef_5': 2004,
        'radius_coef_6': 2005,
        'radius_coef_7': 2006,
        'radius_coef_8': 2007,
    }

    # Map from parameters names to the corresponding key
    _param_map = {
        'glot_pitch': 0,
        'glot_vol': 1,
        'asp_vol': 2,
        'fric_vol': 3,
        'fric_pos': 4,
        'fric_cf': 5,
        'fric_bw': 6,
        'r': (7, _total_sections),
        'r1': 7,
        'r2': 8,
        'r3': 9,
        'r4': 10,
        'r5': 11,
        'r6': 12,
        'r7': 13,
        'r8': 14,
        'velum': 15
    }

    voice_baby = { # XXX this does not work due to a too short vocal tract
        'vocal_tract_length': 7.5,
        'glottal_pulse_tp': 40.0,
        'glottal_pulse_tn_min': 24.0,
        'glottal_pulse_tn_max': 24.0,
        'reference_glottal_pitch': 7.5,
        'breathiness': 1.5,
        'aperture_radius': 3.05,
        'nose_radius': [1.35, 1.96, 1.91, 1.3, 0.73],
        'radius_coef': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    }

    voice_female = {
        'vocal_tract_length': 15.0,
        'glottal_pulse_tp': 40.0,
        'glottal_pulse_tn_min': 32.0,
        'glottal_pulse_tn_max': 32.0,
        'reference_glottal_pitch': 0.0,
        'breathiness': 1.5,
        'aperture_radius': 3.05,
        'nose_radius': [1.35, 1.96, 1.91, 1.3, 0.73],
        'radius_coef': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    }

    voice_large_child = {
        'vocal_tract_length': 12.5,
        'glottal_pulse_tp': 40.0,
        'glottal_pulse_tn_min': 24.0,
        'glottal_pulse_tn_max': 24.0,
        'reference_glottal_pitch': 2.5,
        'breathiness': 1.5,
        'aperture_radius': 3.05,
        'nose_radius': [1.35, 1.96, 1.91, 1.3, 0.73],
        'radius_coef': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    }

    voice_male = {
        'vocal_tract_length': 17.5,
        'glottal_pulse_tp': 40.0,
        'glottal_pulse_tn_min': 24.0,
        'glottal_pulse_tn_max': 24.0,
        'reference_glottal_pitch': -12.0,
        'breathiness': 0.5,
        'aperture_radius': 3.05,
        'nose_radius': [1.35, 1.96, 1.91, 1.3, 0.73],
        'radius_coef': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    }

    voice_small_child = {
        'vocal_tract_length': 10.0,
        'glottal_pulse_tp': 40.0,
        'glottal_pulse_tn_min': 24.0,
        'glottal_pulse_tn_max': 24.0,
        'reference_glottal_pitch': 5.0,
        'breathiness': 0.5,
        'aperture_radius': 3.05,
        'nose_radius': [1.35, 1.96, 1.91, 1.3, 0.73],
        'radius_coef': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    }

    @staticmethod
    def _load_library():
        '''
        Loads the TRM library if this has not already been done.
        '''

        if not TRM._lib is None:
            return  # Library already loaded, all set

        # Iterate over the following set of locations, accept the first
        # working location
        this_dir = os.path.dirname(os.path.realpath(__file__))
        locations = [
            os.path.join(this_dir, '..', 'build', TRM._libname),
            'libgnuspeech_trm.so'
        ]
        for i, location in enumerate(locations):
            try:
                TRM._lib = ctypes.cdll.LoadLibrary(location)
                # Apparently this worked! Continue.
                break
            except OSError:
                # If this is the last try, raise the exception, otherwise, just
                # continue
                if i == len(locations) - 1:
                    raise

        # Set the argument and return types correctly
        TRM._lib.gnuspeech_trm_create.argstype = []
        TRM._lib.gnuspeech_trm_create.restype = TRM._c_void_p

        TRM._lib.gnuspeech_trm_reset.argstype = [TRM._c_void_p]
        TRM._lib.gnuspeech_trm_reset.restype = None

        TRM._lib.gnuspeech_trm_free.argstype = [TRM._c_void_p]
        TRM._lib.gnuspeech_trm_free.restype = None

        TRM._lib.gnuspeech_trm_set_config.argstype = \
               [TRM._c_void_p, ctypes.c_int, ctypes.c_double]
        TRM._lib.gnuspeech_trm_set_config.restype = ctypes.c_bool

        TRM._lib.gnuspeech_trm_get_config.argstype = \
              [TRM._c_void_p, ctypes.c_int, TRM._c_double_p]
        TRM._lib.gnuspeech_trm_get_config.restype = ctypes.c_bool

        TRM._lib.gnuspeech_trm_set_parameter.argstype = \
               [TRM._c_void_p, ctypes.c_int, ctypes.c_double]
        TRM._lib.gnuspeech_trm_set_parameter.restype = ctypes.c_bool

        TRM._lib.gnuspeech_trm_get_parameter.argstype = \
              [TRM._c_void_p, ctypes.c_int, TRM._c_double_p]
        TRM._lib.gnuspeech_trm_get_parameter.restype = ctypes.c_bool

        TRM._lib.gnuspeech_trm_set_parameters.argstype = \
              [TRM._c_void_p, TRM._c_double_p]
        TRM._lib.gnuspeech_trm_set_parameters.restype = ctypes.c_bool

        TRM._lib.gnuspeech_trm_get_parameters.argstype = \
              [TRM._c_void_p, TRM._c_double_p]
        TRM._lib.gnuspeech_trm_get_parameters.restype = ctypes.c_bool

        TRM._lib.gnuspeech_trm_set_parameter_dynamics.argstype = \
              [TRM._c_void_p, TRM._c_double_p]
        TRM._lib.gnuspeech_trm_set_parameter_dynamics.restype = ctypes.c_bool

        TRM._lib.gnuspeech_trm_get_parameter_dynamics.argstype = \
              [TRM._c_void_p, TRM._c_double_p]
        TRM._lib.gnuspeech_trm_get_parameter_dynamics.restype = ctypes.c_bool

        TRM._lib.gnuspeech_trm_synthesize.argstype = \
              [TRM._c_void_p, TRM._c_double_p, ctypes.c_uint, ctypes.c_bool]
        TRM._lib.gnuspeech_trm_synthesize.restype = ctypes.c_int

        return TRM._lib

    def __init__(self, config_dict=None):
        '''
        Allocates a new TRM instance in the C backend.
        '''

        # Make sure the library is loaded
        self._load_library()

        # Create a new instance
        self._inst = TRM._lib.gnuspeech_trm_create()

        # Load configuration from the given dictionary
        if not config_dict is None:
            self.configure(config_dict)

    def __del__(self):
        '''
        Frees the memory allocated in the C backend.
        '''

        TRM._lib.gnuspeech_trm_free(self._inst)

    def __getattr__(self, name):
        '''
        Allows to read a configuration option or parameter as a property of a
        TRM instance.
        '''

        # Make sure the attribute is in the corresponding map
        mapping = None
        if name in TRM._conf_map:
            mapping, f = TRM._conf_map, TRM._lib.gnuspeech_trm_get_config
        elif name in TRM._param_map:
            mapping, f = TRM._param_map, TRM._lib.gnuspeech_trm_get_parameter
        else:
            raise AttributeError()

        # Special handling for array types
        descr = mapping[name]
        if isinstance(descr, tuple):
            offs, n = descr
        else:
            offs, n = descr, 1
        res = [None] * n

        # Try to read the values
        for i in range(n):
            x = ctypes.c_double(0.0)
            ok = f(self._inst, offs + i, ctypes.pointer(x))
            if not ok:
                raise AttributeError()
            res[i] = x.value

        # Unpack single values
        return res[0] if n == 1 else res

    def __setattr__(self, name, value):
        '''
        Allows to set a configuration option or parameter by directly assigning
        it to the TRM instance.
        '''

        # Make sure the attribute is in the corresponding map
        mapping, f = None, None
        if name in TRM._conf_map:
            mapping, f = TRM._conf_map, TRM._lib.gnuspeech_trm_set_config
        elif name in TRM._param_map:
            mapping, f = TRM._param_map, TRM._lib.gnuspeech_trm_set_parameter
        else:
            super(TRM, self).__setattr__(name, value)
            return

        # Special handling for array types
        descr = mapping[name]
        if isinstance(descr, tuple):
            offs, n = descr
            if len(value) != n:
                raise Exception(
                    'value must by an array of length {}'.format(n))
        else:
            offs, n, value = descr, 1, [value]

        # Set all elements
        for i in range(n):
            x = ctypes.c_double(value[i])
            ok = f(self._inst, offs + i, x)
            if not ok:
                raise ValueError()

    def reset(self):
        '''
        Resets this instance to default settings.
        '''
        TRM._lib.gnuspeech_trm_reset(self._inst)

    def configure(self, configuration):
        '''
        Loads configuration options and parameters from the given dictionary.
        '''
        for key, value in configuration.items():
            setattr(self, key, value)

    def set_parameters(self, params):
        if len(params) != TRM._total_params:
            raise Exception('parameters must be a {}-dimensional array'.format(
                TRM._total_params))
        params = np.array(params).flatten().astype(
            np.float64, order='C', copy=False)
        ok = TRM._lib.gnuspeech_trm_set_parameters(
            self._inst, params.ctypes.data_as(TRM._c_double_p))
        if not ok:
            raise Exception('Error while setting parameters')

    def get_parameters(self):
        params = np.empty(TRM._total_params, dtype=np.float64, order='C')
        ok = TRM._lib.gnuspeech_trm_get_parameters(
            self._inst, params.ctypes.data_as(TRM._c_double_p))
        if not ok:
            raise Exception('Error while getting parameters')
        return params

    def set_parameter_dynamics(self, params):
        if len(params) != TRM._total_params:
            raise Exception('parameters must be a {}-dimensional array'.format(
                TRM._total_params))
        params = np.array(params).flatten().astype(
            np.float64, order='C', copy=False)
        ok = TRM._lib.gnuspeech_trm_set_parameter_dynamics(
            self._inst, params.ctypes.data_as(TRM._c_double_p))
        if not ok:
            raise Exception('Error while setting parameters')

    def get_parameter_dynamics(self):
        params = np.empty(TRM._total_params, dtype=np.float64, order='C')
        ok = TRM._lib.gnuspeech_trm_get_parameter_dynamics(
            self._inst, params.ctypes.data_as(TRM._c_double_p))
        if not ok:
            raise Exception('Error while getting parameters')
        return params

    def synthesize(self, n_samples_max=4096, flush=False):
        tar = np.empty(n_samples_max, dtype=np.float32, order='C')
        n = TRM._lib.gnuspeech_trm_synthesize(
            self._inst, tar.ctypes.data_as(TRM._c_float_p), n_samples_max,
            flush)
        if n < 0:
            raise Exception('Cannot initialize synthesizer')
        return tar[:n]


if __name__ == '__main__':
    import sys
    trm = TRM(TRM.voice_female)

    # Test reading/writing configuration and parameter entries
    for key in trm._conf_map.keys():
        value = trm.__getattr__(key)
        print(key, value)
        trm.__setattr__(key, value)
    for key in trm._param_map.keys():
        value = trm.__getattr__(key)
        trm.__setattr__(key, value)

    # Generate a standard, relaxed open mouth tone ('äääääh')
    trm.volume = 60.0
    trm.filter_period = 10e-3
    trm.glot_pitch = 0.0
    trm.glot_vol = 60.0
    trm.velum = 0.0
    trm.fric_vol = 0.0
    trm.fric_pos = 0.0
    trm.fric_bw = 220.0
    trm.fric_cf = 100.0
    trm.r = [1.0] * trm._total_sections
    data = trm.synthesize(int(trm.output_rate * 2))
    sys.stdout.buffer.write(data.tobytes())

