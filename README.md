# LoRaFREE

This is repository of LoRaFREE, a new simulator tool that have been used in [1]. Please cite the paper if you used LoRaFREE in your research.

LoRaFREE is a more comprehensive Simpy simulator for LoRa than LoRaSim, which considers a packet error model, the imperfect orthogonality of spreading factors, the fading impact, and the duty cycle limitation at both, the devices and the gateway. Moreover, LoRaFREE supports bidirectional communication by adding the downlink capability and the retransmission strategy in case of confirmable uplink transmissions. Furthermore, we extended the energy consumption profile from LoRaSim to consider the consumed energy at the reception time. In additions to that, you will find a simulator for the synchronized transmissions schedule based on the work in [1].


[1] @ARTICLE{abdelfadeelfree2018,
author={K. Q. {Abdelfadeel} and D. {Zorbas} and V. {Cionca} and D. {Pesch}},
journal={IEEE Internet of Things Journal},
title={FREE -Fine-grained Scheduling for Reliable and Energy Efficient Data Collection in LoRaWAN},
year={2019},
volume={},
number={},
pages={1-1},
keywords={Logic gates;Synchronization;Data collection;Schedules;Reliability;Uplink;Performance evaluation},
doi={10.1109/JIOT.2019.2949918},
ISSN={2372-2541},
month={},
}
