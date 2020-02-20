import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PastaPhase", # Replace with your own username
    version="0.0.1",
    author="P. P. Galuzio and G. Grams",
    author_email="galuzio.paulo@protonmail.com",
    description="Neutron Star package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ppgaluzio/PastaPhase",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
