# MGcount
MGcount is a programm for counting whole-transcriptome RNA-seq reads. 

## Requirements
MGcount deppends on FeatureCounts. Please download it from the following link: http://subread.sourceforge.net/

## Installation
MGcount is written in Python and is executed from the command line. You can either download the executable version (single binary file) or install it as a Python3 module. 

### Install and run as a Python3 module
You can install the package as a Python module:

```shell
pip3 install git+https://github.com/hitaandrea/MGcount.git
```
Once the package is installed, run the tool as a Python installed module:

```shell
python3 -m mgcount [args]
```

### Download and run the executable program
Alternatively can download the latest release as a binary executable file [here](https://github.com/hitaandrea/MGcount/releases/download/1.0.0-beta/MGcount). 

Save the program file to your Linux system and set the permissions to allow executing the file as program:

```shell
chmod +x mgcount
```

Once the file is executable, run the tool by calling the file from the command line with the desired arguments.

```shell
MGcount [args]
```

## Documentation
MGcount user guide can be found inside the docs subfolder or accessed through the following link:
[Link to User Guide!](https://filedn.com/lTnUWxFTA93JTyX3Hvbdn2h/mgcount/UserGuide.html)


## Credits
The work is funded by a Marie Curie early stage researcher fellowship. (European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No. 813282).

## Licences
MGcount itself is free software distributed under GPL.
