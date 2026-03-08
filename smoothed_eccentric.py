import numpy as np
from scipy.interpolate import PchipInterpolator

def segment_average_broadcast(data, indices):
    """
    Replaces values in 'data' with the average of their segment, 
    maintaining the original shape.
    """
    # 1. Ensure we cover the full range of the array
    # We want breakpoints at 0, your indices, and the end of the array
    breakpoints = np.unique(np.concatenate(([0], indices, [len(data)])))
    
    # 2. Sum each segment (O(N))
    # Note: reduceat takes the start index of each segment
    sums = np.add.reduceat(data, breakpoints[:-1])
    
    # 3. Calculate lengths of each segment
    counts = np.diff(breakpoints)
    
    # 4. Calculate averages
    averages = sums / counts
    
    # 5. Expand (Repeat) averages to the foriginal shape (O(N))
    return averages, np.repeat(averages, counts)

def fit_piecewise_monotonic(x, y):
    """
    Fits a piecewise cubic polynomial that is guaranteed 
    to preserve monotonicity (no ripples/overshoot).
    """
    # Sort data by x to ensure the interpolator works correctly
    idx = np.argsort(x)
    x_sorted = x[idx]
    y_sorted = y[idx]
    
    # PchipInterpolator creates the piecewise function
    # It calculates derivatives at knots to ensure smoothness 
    # and monotonicity.
    pcs = PchipInterpolator(x_sorted, y_sorted)
    
    return pcs


def smoothed_waveform(hp, hc):
    import numpy
    amp = (hp**2.0 + hc**2.0) ** 0.5
    amp2 = amp.copy()
    
    phase = amp.copy()
    phase.data[:] = numpy.unwrap(numpy.arctan2(hp, hc))
    phase -= phase[0]
    phase2 = phase.copy()
    
    # Average the phase evolution over one orbital cycle / 2 GW cycles
    ncycle = 2

    mmax = numpy.argmax(amp)
    nfull = int((len(amp) - mmax) * 2)
    nfull2 = int(len(amp) - numpy.where(amp > amp[mmax] * 0.8)[0][0])
    nfull = max(nfull, nfull2)
    
    #nfull = 300 
    nstart = 0

    valid = slice(nstart, -nfull)
    rval = numpy.arange(-phase[nstart], -phase[-nfull], numpy.pi * 2 * ncycle)
    points = numpy.searchsorted(-phase[valid], rval)
    points[-1] -= 1
    
    ptimes = phase.sample_times[valid]
    
    avg, a2 = segment_average_broadcast(amp[valid], points)
    
    sa = numpy.interp(ptimes, ptimes[points[:len(avg)]], avg)
    
    amp2.data[valid] = sa
    phase2.data[valid] = fit_piecewise_monotonic(ptimes[points], phase[valid][points])(ptimes)
    
    # Reconstruct new smoothed template
    h =  (numpy.exp(1.0j * (phase2 - np.pi)) * amp2)
    return h.real(), h.imag()

def dominant_harmonic_waveform(**params):
    from pycbc.waveform import get_td_waveform
    if 'approximant' in params:
        params.pop('approximant')

    base = params['base_model']
    if hasattr(base, 'decode'):
        base = base.decode()

    hp, hc = get_td_waveform(approximant=base, **params)
    return smoothed_waveform(hp, hc)
