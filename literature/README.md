# Literature Review

This directory contains the literature that shaped my understanding of quantum key distribution (QKD), Twin-Field QKD (TF-QKD), measurement-device-independent protocols, satellite quantum communication, and the theoretical limits of quantum communications.

The goal of this collection is not simply to archive papers, but to record the concepts that influenced my research.

---

## Papers Reviewed

### 1. Overcoming the Rate–Distance Limit of Quantum Key Distribution without Quantum Repeaters
**Authors:** Marco Lucamarini et al. (Nature, 2018)

**Focus**
- Original Twin-Field QKD proposal.
- Introduces the twin-field concept.
- Shows key-rate scaling proportional to √η instead of η.
- Uses an untrusted central measurement station.

**Takeaways**
- Learned why TF-QKD can surpass repeaterless rate-distance limits.
- Understood the importance of single-photon interference.
- Introduced me to phase slicing and twin-field matching.

**Relation to my work**
Foundation paper for understanding TF-QKD.

---

### 2. Sending or Not Sending: Twin-Field Quantum Key Distribution
**Authors:** Xiang-Bin Wang et al.

**Focus**
- SNS-TF-QKD protocol.
- Addresses security issues present in the original TF-QKD proposal.
- Improves tolerance to misalignment errors.

**Takeaways**
- Learned how the "sending or not sending" idea simplifies security analysis.
- Understood why the original decoy-state formulas cannot always be directly applied.
- Saw how SNS improves practical performance over long distances.

**Relation to my work**
Reference for practical TF-QKD implementations.

---

### 3. Phase-Matching Quantum Key Distribution (PM-QKD)
**Authors:** Xiongfeng Ma et al.

**Focus**
- Optical-mode security proof.
- Measurement-device-independent protocol.
- Phase post-compensation instead of phase locking.

**Takeaways**
- Learned the optical-mode security framework.
- Understood phase post-compensation techniques.
- Saw another protocol achieving √η scaling.

**Relation to my work**
Provides theoretical insight into TF-QKD security.

---

### 4. Making the Decoy-State MDI-QKD Practically Useful
**Authors:** Yi-Heng Zhou et al.

**Focus**
- Practical improvements for decoy-state MDI-QKD.
- Four-intensity protocol.
- Better finite-size performance.

**Takeaways**
- Learned practical parameter optimization.
- Better understanding of decoy-state estimation.
- Importance of statistical fluctuation analysis.

**Relation to my work**
Useful background before studying TF-QKD.

---

### 5. Micius Quantum Experiments in Space
**Authors:** Chao-Yang Lu et al.

**Focus**
- Review of satellite quantum communication experiments.
- Micius satellite missions.
- Space-ground QKD.
- Entanglement distribution.
- Quantum teleportation.

**Takeaways**
- Learned how satellite QKD has evolved experimentally.
- Understood advantages of satellite channels over long optical fibers.
- Gained insight into practical engineering challenges.

**Relation to my work**
Reference for satellite-based quantum communication.

---

### 6. PLOB Bound in Fiber vs Free Space
**Author:** Stefano Pirandola (Seminar Slides)

**Focus**
- Fundamental limits of quantum communication.
- PLOB bound.
- Fiber vs free-space communication.
- Repeaterless quantum communication.

**Takeaways**
- Learned why the PLOB bound is the benchmark for repeaterless communication.
- Understood how free-space channels differ from fiber channels.
- Better interpretation of theoretical communication limits.

**Relation to my work**
Used as the theoretical benchmark when evaluating QKD protocols.

---

### 7. Implementation Attacks against QKD Systems
**Published by:** German Federal Office for Information Security (BSI)

**Focus**
- Practical attacks against QKD implementations.
- Security assumptions.
- Device vulnerabilities.
- Countermeasures.

**Takeaways**
- Learned that practical implementations introduce vulnerabilities beyond protocol security.
- Better understanding of implementation security.
- Importance of considering real devices in secure system design.

**Relation to my work**
Reminder that protocol security and implementation security are different problems.

---

## Overall Learning

Reading these papers helped me understand

- Evolution from BB84 to MDI-QKD and TF-QKD.
- Why repeaterless communication has theoretical limits.
- How TF-QKD changes key-rate scaling.
- Practical protocol improvements for long-distance QKD.
- Satellite-based quantum communication.
- Security considerations beyond theoretical proofs.

---

## Research Areas Covered

- Quantum Key Distribution
- Twin-Field QKD
- SNS-TF-QKD
- PM-QKD
- MDI-QKD
- Decoy-State Methods
- Satellite Quantum Communication
- PLOB Bound
- Quantum Communication Security