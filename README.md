Project Aim:
The primary objective of this project is to design and implement a Quadrature Phase Shift Keying (QPSK) transceiver on an FPGA platform. The project will be executed in two main phases, each focusing on key aspects of digital communication system design and validation.

Phase 1: MATLAB-Based Logic Simulation and System Design
In the initial phase, MATLAB will be utilized for the high-level modeling and simulation of the QPSK transceiver. This includes:

Designing and optimizing transmitter and receiver filters (e.g., Root Raised Cosine filters).

Implementing and verifying QPSK symbol mapping and demapping logic.

Simulating baseband signal processing to validate theoretical performance before hardware implementation.

This phase serves as a functional proof-of-concept and allows for precise tuning of signal processing algorithms before transitioning to the hardware domain.

Phase 2: RTL Development and FPGA Implementation
The second phase focuses on RTL-level design and implementation of the system on an FPGA. This includes:

Developing the complete digital signal chain for both transmitter and receiver.

Implementing the CORDIC algorithm for efficient carrier phase generation and demodulation.

Designing and integrating preamble detection logic for robust frame synchronization.

Implementing time synchronization to align received symbols accurately in time.

Incorporating channel emulation to test the robustness of the transceiver under various noise and interference conditions.

The final outcome will be a functional, real-time QPSK transceiver system implemented on FPGA, demonstrating both the theoretical and practical aspects of modern digital communication design.
![QPSK_Transiver](https://github.com/user-attachments/assets/d2dea7f3-093a-436b-8464-0b25a0ec3a1b)
