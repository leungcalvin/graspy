import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="graspy", # Replace with your own username
    version="0.0.1",
    author="Calvin Leung",
    author_email="calvinl@mit.edu",
    description="A Python interface to GRASP 2018",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/leungcalvin/graspy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
