# LoRaFREE

This is repository of LoRaFREE, a new simulator tool that have been used in [1]. Please cite the paper if you consider using LoRaFREE in your research.

LoRaFREE is a more comprehensive Simpy simulator for LoRa than LoRaSim, which considers a packet error model, the imperfect orthogonality of spreading factors, the fading impact, and the duty cycle limitation at both, the devices and the gateway. LoRaFREE supports bidirectional communication by adding the downlink capability and the retransmission strategy in case of confirmable uplink transmissions. LoRaFree also extends the energy consumption profile from LoRaSim to consider the consumed energy at the reception time. In additions to that, we added an additional simulator for the synchronized transmissions schedule based on the work in [1].


[1] @ARTICLE{8884111,
  author={K. Q. {Abdelfadeel} and D. {Zorbas} and V. {Cionca} and D. {Pesch}},
  journal={IEEE Internet of Things Journal}, 
  title={ $FREE$ —Fine-Grained Scheduling for Reliable and Energy-Efficient Data Collection in LoRaWAN}, 
  year={2020},
  volume={7},
  number={1},
  pages={669-683},}
