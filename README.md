# MultiLoc

MultiLoc is a MutiBand Localization system that merges information from mmWave and sub-6GHz bands. MultiLoc is implemented on Commercial-off-the-shelf  (COST) devices and is the first system to achieve less than 17 cm median location error with (COST) under their normal mode of operation.

## Table of content

- [Getting started](#getting-started)
- [MultiLoc implementation](#multiloc-implementation)

# Getting started

MultiLoc is implemented in 2 sub-6GHz devices (CSI and FTM) and one mmWave device (Mikrotik) as depicted in the image below:

<img src="https://github.com/IMDEANetworksWNG/MultiLoc/blob/main/implementation.png" width="400" height="400">

You can find the intructions below to get the CSI and FTM extractor tools for all the devices below.

**Sub-6GHZ devices**

* CSI  router (ASUS AC1300)

The intructions for running CSI extractor tool are found in this [link](https://github.com/IMDEANetworksWNG/UbiLocate)

* FTM  router (ASUS AC1300)

The intructions for running FTM are found in this [link](https://www.winlab.rutgers.edu/~gruteser/projects/ftm/index.htm)


**mmWave device**

* Mikrotik 60GHz device

The intructions for running CSI and FTM extractor tools are found in this [link](https://github.com/IMDEANetworksWNG/Mikrotik-researcher-tools)

In addition, this [link](https://github.com/IMDEANetworksWNG/Mikrotik-researcher-tools) contains the intruction for calibrating the Mikrotik devices. 

# MultiLoc implementation

