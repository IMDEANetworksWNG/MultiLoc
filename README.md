# MultiLoc

MultiLoc is a MutiBand Localization system that merges information from mmWave and sub-6GHz bands. MultiLoc is implemented on Commercial-Off-The-Shelf  (COST) devices and is the first system to achieve less than 18 cm median location error with COST under their normal mode of operation.

If you find MultiLoc useful, we kindly ask you to cite our paper:
```
@inproceedings{MultiLoc,
author = {Blanco, Alejandro and Mateo, Pablo Jim\'{e}nez and Gringoli, Francesco and Widmer, Joerg},
title = {Augmenting MmWave Localization Accuracy through Sub-6 GHz on off-the-Shelf Devices},
year = {2022},
isbn = {9781450391856},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3498361.3538920},
doi = {10.1145/3498361.3538920},
booktitle = {Proceedings of the 20th Annual International Conference on Mobile Systems, Applications and Services},
pages = {477–490},
numpages = {14},
keywords = {AoA, indoor localization, wireless networks, ToF, mmWave, CSI},
location = {Portland, Oregon},
series = {MobiSys '22}
}
```

## Table of content

- [Getting started](#getting-started)
- [MultiLoc implementation](#multiloc-implementation)

# Getting started

MultiLoc is implemented in 2 sub-6GHz devices (one for CSI and another for FTM) and in the MikroTik 60GHz device as depicted in the image below:

<img src="https://github.com/IMDEANetworksWNG/MultiLoc/blob/main/implementation.png" width="400" height="400">

You can find the instructions below to get the CSI and FTM extractor tools for all the devices.

**Sub-6GHZ devices**

* CSI  router (ASUS AC1300). The instructions for running the CSI extractor tool are found in this [link](https://github.com/IMDEANetworksWNG/UbiLocate).

* FTM  router (ASUS AC2900). The instructions for running FTM are found in this [link](https://www.winlab.rutgers.edu/~gruteser/projects/ftm/index.htm).


**mmWave device**

* Mikrotik 60GHz device. The instructions for running CSI and FTM extractor tools are found in this [link](https://github.com/IMDEANetworksWNG/Mikrotik-researcher-tools). In addition, this [link](https://github.com/IMDEANetworksWNG/MikroTik-mD-Track) contains the instructions for calibrating the Mikrotik devices. 

# MultiLoc implementation
Run matlab_scripts/Main.m to run the main of the MultiLoc implementation.

Take into account that the script Main_HF.m is commented since it cleans the mmWave CSI. The CSI data cleaning uses GMM and we believe that every time is executed, it does not produce exactly same output. We recommend using the already processed mmWave data of the project.

The scripts Main_LFs are also commented since the mD-track algorithm takes a considerable amount of time to get the path parameters, but feel free to execute it.
