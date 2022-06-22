# MultiLoc

MultiLoc is a MutiBand Localization system that merges information from mmWave and sub-6GHz bands. MultiLoc is implemented on Commercial-Off-The-Shelf  (COST) devices and is the first system to achieve less than 17 cm median location error with COST under their normal mode of operation.

## Table of content

- [Getting started](#getting-started)
- [MultiLoc implementation](#multiloc-implementation)

# Getting started

MultiLoc is implemented in 2 sub-6GHz devices (CSI and FTM) and one mmWave device (Mikrotik) as depicted in the image below:

<img src="https://github.com/IMDEANetworksWNG/MultiLoc/blob/main/implementation.png" width="400" height="400">

You can find the instructions below to get the CSI and FTM extractor tools for all the devices below.

**Sub-6GHZ devices**

* CSI  router (ASUS AC1300). The instructions for running the CSI extractor tool are found in this [link](https://github.com/IMDEANetworksWNG/UbiLocate).

* FTM  router (ASUS AC2900). The instructions for running FTM are found in this [link](https://www.winlab.rutgers.edu/~gruteser/projects/ftm/index.htm).


**mmWave device**

* Mikrotik 60GHz device. The instructions for running CSI and FTM extractor tools are found in this [link](https://github.com/IMDEANetworksWNG/Mikrotik-researcher-tools). In addition, this [link](https://github.com/IMDEANetworksWNG/MikroTik-mD-Track) contains the instructions for calibrating the Mikrotik devices. 

# MultiLoc implementation (dataset + localization algorithms)
Run matlab_scripts/Main.m to run the main of the MultiLoc implementation.

Take into account that the Main_HF.m is commented since it extracts the mmWave CSI and ToF. The CSI data cleaning uses GMM and we believe that every time is executed it does not produce the same output. We recomment yo use the already processed data of the project.

The Main_LFs are also commented since mD-track takes time to get the path parameters, but feel free to execute it.
