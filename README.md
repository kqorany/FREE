# LoRaFREE

This repository contains LoRaFREE, a new simulator tool that have been used in [1]. Please cite the paper if you used LoRaFREE in your research.

LoRaFREE is a more comprehensive Simpy simulator for LoRa than LoRaSim, which considers a packet error model, the imperfect orthogonality of spreading factors, the fading impact, and the duty cycle limitation at both, the devices and the gateway. Moreover, LoRaFREE supports bidirectional communication by adding the downlink capability and the retransmission strategy in case of confirmable uplink transmissions. Furthermore, we extended the energy consumption profile from LoRaSim to consider the consumed energy at the reception time. In additions to that, you will find a simulator for the synchronized transmissions schedule based on the work in [1].


[1] @article{abdelfadeel2018free,
  title={FREE-Fine-grained Scheduling for Reliable and Energy Efficient Data Collection in LoRaWAN},
  author={Abdelfadeel, Khaled Q and Zorbas, Dimitrios and Cionca, Victor and O'Flynn, Brendan and Pesch, Dirk},
  journal={arXiv preprint arXiv:1812.05744},
  year={2018}
}
