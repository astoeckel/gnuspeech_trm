/**
 *  gnuspeech_trm -- Standalone gnuspeech articulatory synthesiser
 *  Copyright (C) 2018  Andreas Stöckel
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef _GNUSPEECH_TRM_H_
#define _GNUSPEECH_TRM_H_

#include <stdbool.h>

// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
#define TRM_DLL_EXPORT __declspec(dllexport)
#else
#if __GNUC__ >= 4
#define TRM_DLL_EXPORT __attribute__((visibility("default")))
#else
#define TRM_DLL_EXPORT
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define TRM_CONF_OUTPUT_RATE 10
#define TRM_CONF_CONTROL_RATE 20
#define TRM_CONF_FILTER_PERIOD 30
#define TRM_CONF_VOLUME 40
#define TRM_CONF_WAVEFORM 50
#define TRM_CONF_TEMPERATURE 60
#define TRM_CONF_LOSS_FACTOR 70
#define TRM_CONF_MOUTH_COEF 80
#define TRM_CONF_NOSE_COEF 90
#define TRM_CONF_THROAT_CUTOFF 100
#define TRM_CONF_THROAT_VOL 110
#define TRM_CONF_MODULATION 120
#define TRM_CONF_MIX_OFFSET 130
#define TRM_CONF_GLOTTAL_PULSE_TP 140
#define TRM_CONF_GLOTTAL_PULSE_TN_MIN 150
#define TRM_CONF_GLOTTAL_PULSE_TN_MAX 160
#define TRM_CONF_BREATHINESS 170
#define TRM_CONF_VOCAL_TRACT_LENGTH 180
#define TRM_CONF_REFERENCE_GLOTTAL_PITCH 190
#define TRM_CONF_APERTURE_RADIUS 200

#define TRM_CONF_NOSE_RADIUS_1 1000
#define TRM_CONF_NOSE_RADIUS_2 1001
#define TRM_CONF_NOSE_RADIUS_3 1002
#define TRM_CONF_NOSE_RADIUS_4 1003
#define TRM_CONF_NOSE_RADIUS_5 1004
#define TRM_CONF_NOSE_RADIUS_6 1005

#define TRM_CONF_RADIUS_COEF_1 2000
#define TRM_CONF_RADIUS_COEF_2 2001
#define TRM_CONF_RADIUS_COEF_3 2002
#define TRM_CONF_RADIUS_COEF_4 2003
#define TRM_CONF_RADIUS_COEF_5 2004
#define TRM_CONF_RADIUS_COEF_6 2005
#define TRM_CONF_RADIUS_COEF_7 2006
#define TRM_CONF_RADIUS_COEF_8 2007

#define TRM_CONF_NOSE_RADIUS_MIN TRM_CONF_NOSE_RADIUS_1
#define TRM_CONF_NOSE_RADIUS_MAX TRM_CONF_NOSE_RADIUS_6
#define TRM_CONF_NOSE_RADIUS_COUNT \
	(TRM_CONF_NOSE_RADIUS_MAX - TRM_CONF_NOSE_RADIUS_MIN + 1)
#define TRM_CONF_RADIUS_COEF_MIN TRM_CONF_RADIUS_COEF_1
#define TRM_CONF_RADIUS_COEF_MAX TRM_CONF_RADIUS_COEF_8
#define TRM_CONF_RADIUS_COEF_COUNT \
	(TRM_CONF_RADIUS_COEF_MAX - TRM_CONF_RADIUS_COEF_MIN + 1)

#define TRM_PARAM_GLOT_PITCH 0
#define TRM_PARAM_GLOT_VOL 1
#define TRM_PARAM_ASP_VOL 2
#define TRM_PARAM_FRIC_VOL 3
#define TRM_PARAM_FRIC_POS 4
#define TRM_PARAM_FRIC_CF 5
#define TRM_PARAM_FRIC_BW 6
#define TRM_PARAM_R1 7
#define TRM_PARAM_R2 8
#define TRM_PARAM_R3 9
#define TRM_PARAM_R4 10
#define TRM_PARAM_R5 11
#define TRM_PARAM_R6 12
#define TRM_PARAM_R7 13
#define TRM_PARAM_R8 14
#define TRM_PARAM_VELUM 15

#define TRM_PARAM_COUNT 16

/**
 * Opaque pointer at a gnuspeech_trm instance.
 */
typedef void *gnuspeech_trm_t;

/**
 * Creates a new gnuspeech_trm instance and returns a pointer at that instance.
 * If memory allocation fails, returns NULL.
 */
TRM_DLL_EXPORT gnuspeech_trm_t gnuspeech_trm_create();

/**
 * Destroys an instance of gnuspeech_trm created previously by
 * gnuspeech_trm_create().
 */
TRM_DLL_EXPORT void gnuspeech_trm_free(gnuspeech_trm_t inst);

/**
 * Resets the gnuspeech instance to its initial state.
 */
TRM_DLL_EXPORT void gnuspeech_trm_reset(gnuspeech_trm_t inst);

/**
 * Sets a synthesiser configuration option. Returns true if the operation was
 * successful, false if either the parameter index or the value was out of
 * bounds.
 */
TRM_DLL_EXPORT bool gnuspeech_trm_set_config(gnuspeech_trm_t inst, int conf,
                                             double value);

/**
 * Gets a synthesiser configuration option and writes it to the given variable.
 * Returns true if the operation was successful, false if the parameter index
 * was out of bounds.
 */
TRM_DLL_EXPORT bool gnuspeech_trm_get_config(const gnuspeech_trm_t inst,
                                             int conf, double *value);

/**
 * Sets a synthesizer parameter to the given value. Returns true if the
 * operation was successful.
 */
TRM_DLL_EXPORT bool gnuspeech_trm_set_parameter(gnuspeech_trm_t inst, int param,
                                                double value);

/**
 * Returns the current value of the given parameter by writing it to the given
 * variable. Returns ture if hte operation was successful.
 */
TRM_DLL_EXPORT bool gnuspeech_trm_get_parameter(const gnuspeech_trm_t inst,
                                                int param, double *value);

/**
 * Reads all parameter values from the given array.
 */
TRM_DLL_EXPORT bool gnuspeech_trm_set_parameters(gnuspeech_trm_t inst,
                                                const double *dx);

/**
 * Writes the current parameter values to the given array.
 */
TRM_DLL_EXPORT bool gnuspeech_trm_get_parameters(gnuspeech_trm_t inst,
                                                 double *dx);

/**
 * Updates the dynamics for all parameters, i.e. the rate with which the
 * parameters should change while synthesizing.
 */
TRM_DLL_EXPORT bool gnuspeech_trm_set_parameter_dynamics(gnuspeech_trm_t inst,
                                                         const double *dx);

/**
 * Writes the dynamics for all parameters to the given array.
 */
TRM_DLL_EXPORT bool gnuspeech_trm_get_parameter_dynamics(gnuspeech_trm_t inst,
                                                         double *dx);

/**
 * Synthesizes the given number of samples and writes them to the given output
 * buffer. Returns the number of samples that actually were written.
 */
TRM_DLL_EXPORT int gnuspeech_trm_synthesize(gnuspeech_trm_t inst,
                                            float *sample_buf,
                                            unsigned int n_samples, bool flush);

#ifdef __cplusplus
}
#endif

#endif /* _GNUSPEECH_TRM_H_ */
