import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="treemer-takkoona",
    version="0.0.1",
    author="Chengyang Ji",
    author_email="chengyang.ji12@alumni.xjtlu.edu.cn",
    description="Phylogeny-dependent method for trimming sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Takkoona/treemer",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
