import setuptools

setuptools.setup(
    name="PaSTA",
    version="0.1",
    author="Ema Smith",
    author_email="ema.smith@yale.edu",
    description="A Python package for plotting time curves and extracting P/S - wave residual ratio. Name stands for P And S wave Time Analysis",
    packages=setuptools.find_packages(include=['PaSTA','PaSTA.*']),
    python_requires='>=3',
    install_requires=["numpy","matplotlib","pandas","obspy","scipy"] 
)
