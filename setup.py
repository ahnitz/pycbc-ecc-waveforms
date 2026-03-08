from setuptools import setup

setup(
    name="pycbc-ecc-waveforms",
    version="0.1.0",
    description="Custom PyCBC plugin for smoothed eccentric waveforms",
    author="Alexander Harevey Nitz",
    py_modules=["smoothed_eccentric"],
    install_requires=[
        "pycbc",
        "numpy",
        "scipy"
    ],
    entry_points={
        'pycbc.waveform.td': [
            # 'name_used_in_flags = module_name:function_name'
            'smoothed_eccentric = smoothed_eccentric:dominant_harmonic_waveform',
        ],
    },
)
