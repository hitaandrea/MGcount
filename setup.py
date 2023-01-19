import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name="mgcount",
    version="1.1.0",
    description="RNA-seq counting tool",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/hitaandrea/MGcount",
    author="Andrea Hita",
    author_email="hitaandrea@gmail.com",
    license="GNU GENERAL PUBLIC LICENSE",
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
    ],
    packages=["mgcount"],
    include_package_data=True,
    package_data={"mgcount": ["data/*.csv"]},
    install_requires=["pandas", "numpy", "scipy", "python-igraph", "pysam", "gtfparse==1.2.1"],
    entry_points={"console_scripts": ["mgcount=mgcount.__main__:main"]},
)

