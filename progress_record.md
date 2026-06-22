# DRDO-SAG Summer Internship

## Project: Analysis of the Performance of Quantum Key Distribution Systems over Free Space Optical Channels and Security Assessment under Quantum Attacks

**Name:** Anand Kumar
**Supervisor:** Mr. Rajesh Kumar, Scientist 'E'
**Organization:** Scientific Analysis Group (SAG), DRDO

---

# Week 01

**Reporting Period:** 16 June 2026 – 21 June 2026

## Objectives

* Understand project objectives and expected deliverables.
* Review existing SNS-TFQKD Optical System Design (OSD).
* Study practical implementations of SNS-TFQKD protocols.
* Develop a simulation framework for protocol performance analysis.
* Investigate Free Space Optical (FSO) channel integration.

---

## 16 June 2026

### Activities Performed

* Attended introductory meeting with the supervisor.
* Discussed internship objectives, expected outcomes, and project roadmap.
* Presented independent work completed prior to internship commencement.
* Received the existing SNS-TFQKD Optical System Design (OSD) for analysis.
* Discussed possible integration of realistic FSO channel effects into the current system.
* Identified major performance parameters relevant to practical QKD deployment.

### Outcome

* Obtained project requirements and simulation resources.
* Established initial understanding of the provided system architecture.

---

## 17 June 2026

### Activities Performed

* Conducted detailed examination of the SNS-TFQKD Optical System Design.
* Reviewed transmitter, channel, and receiver architectures.
* Investigated the implementation methodology adopted in the existing model.
* Studied practical challenges associated with realistic quantum communication systems.
* Reviewed literature related to satellite quantum communication and free-space quantum links.

### Outcome

* Acquired familiarity with the provided simulation architecture.
* Identified key areas requiring further investigation.

---

## 18 June 2026

### Activities Performed

* Studied the paper:
  *"Sending-or-Not-Sending Twin-Field Quantum Key Distribution in Practice"* by Zong-Wen Yu et al.
* Examined practical SNS-TFQKD implementation procedures.
* Studied decoy-state methodology and security analysis.
* Investigated gain estimation, yield estimation, and secure key rate calculations.
* Analyzed detector imperfections, channel losses, and misalignment effects.
* Initiated implementation of protocol equations using Python.

### Outcome

* Developed understanding of practical SNS-TFQKD operation.
* Established theoretical foundation for protocol simulation.

---

## 19 June 2026

### Activities Performed

* Continued implementation of SNS-TFQKD protocol equations.
* Developed channel transmittance models.
* Implemented detector efficiency and dark count models.
* Studied yield and gain calculations.
* Performed preliminary testing of protocol simulation modules.

### Outcome

* Generated initial simulation framework for protocol analysis.

---

## 20 June 2026

### Activities Performed

* Extended simulation modules for protocol performance evaluation.
* Conducted preliminary analysis of secure key rate behavior.
* Investigated the effect of channel loss and misalignment errors.
* Compared simulation outputs with theoretical expectations.

### Outcome

* Generated initial protocol performance results.
* Established baseline framework for future FSO integration.

---

## Literature Survey Conducted During Week 01

* BSI Report: *QKD Systems*
* *Micius Quantum Experiments in Space*
* SNS-TFQKD literature by Zong-Wen Yu et al.
* Twin-Field QKD fundamentals.
* Decoy-state methods.
* Secure key rate estimation techniques.
* Satellite-based quantum communication systems.

---

## Progress Achieved

* Completed preliminary analysis of SNS-TFQKD architecture.
* Developed initial Python simulation framework.
* Implemented detector and channel models.
* Generated preliminary performance results.
* Identified requirements for realistic FSO integration.

---

## Challenges Encountered

* Understanding security proofs and decoy-state analysis.
* Translating protocol equations into simulation models.
* Validation of theoretical assumptions against implementation requirements.

---

## Plan for Subsequent Work

* Validate simulation results with published literature.
* Incorporate realistic FSO channel effects.
* Study atmospheric attenuation and turbulence.
* Continue SNS-TFQKD performance analysis.

---

# Week 02

**Reporting Period:** 22 June 2026 – 28 June 2026

## Objectives

* Analyze the provided SNS-TFQKD OptiSystem implementation.
* Understand OptiSystem architecture and component interactions.
* Study MATLAB-OptiSystem integration.
* Identify protocol-level implementation gaps.
* Develop Send/Not-Send functionality within the OptiSystem framework.

---

## 22 June 2026

### Activities Performed

* Studied OptiSystem documentation, user manuals, and component libraries.
* Examined OptiSystem subsystem architecture and signal processing workflow.
* Investigated MATLAB component integration and co-simulation methodology.
* Reviewed the provided SNS-TFQKD OptiSystem design.

#### Alice Station Analysis

* Examined CW laser source configuration.
* Studied modulation and encoding blocks.
* Investigated phase randomization implementation.
* Reviewed intensity modulation and optical transmission chain.

#### Bob Station Analysis

* Reviewed transmitter architecture.
* Compared implementation with Alice station.
* Examined modulation and channel configuration.

#### Charlie Station Analysis

* Studied interference setup using optical coupler.
* Investigated single-photon detector arrangement.
* Reviewed detector output processing blocks.

### Protocol Assessment

* Compared the provided design against standard SNS-TFQKD protocol requirements.
* Identified reliance on User Defined Bit Sequence Generators.
* Investigated limitations of deterministic state preparation.
* Examined protocol compliance regarding Send/Not-Send operation.
* Investigated absence of explicit vacuum-state generation.
* Reviewed implementation requirements for protocol-level randomness.

### MATLAB Integration Work

* Studied custom MATLAB component implementation within OptiSystem.
* Investigated optical and electrical signal interfaces.
* Began development of a MATLAB-based Send/Not-Send decision module.
* Explored methods for generating vacuum states corresponding to NOT-SEND events.
* Evaluated pseudo-random signal generation for protocol control.

### Literature and Documentation Studied

* SNS-TFQKD protocol fundamentals.
* Twin-Field QKD architecture.
* Measurement-Device-Independent QKD concepts.
* OptiSystem user documentation.
* MATLAB-OptiSystem co-simulation documentation.
* Send/Not-Send state preparation methodologies.
* Vacuum-state encoding concepts.

### Outcome

* Completed initial review of the provided OptiSystem architecture.
* Identified protocol-level implementation gaps.
* Acquired working familiarity with OptiSystem workflow.
* Established implementation strategy for introducing Send/Not-Send functionality.
* Prepared foundation for protocol-compliant modification of the existing model.

---

## Progress Achieved (Up to 22 June 2026)

* Completed review of the provided SNS-TFQKD OptiSystem design.
* Studied OptiSystem architecture and MATLAB integration workflow.
* Identified protocol-level deficiencies in the existing implementation.
* Initiated development of a Send/Not-Send decision mechanism.
* Established a roadmap for further protocol enhancement.

---

## Challenges Encountered

* Limited visibility into internal subsystem implementations.
* Understanding OptiSystem-specific signal structures.
* Mapping theoretical SNS-TFQKD operations into OptiSystem components.
* Verification of protocol compliance within the provided design.

---

## Planned Activities

* Complete MATLAB-based Send/Not-Send module implementation.
* Validate integration with OptiSystem.
* Analyze phase randomization subsystem in detail.
* Investigate decoy-state implementation.
* Develop detection statistics extraction methodology.
* Implement gain, QBER, and secret key rate calculations.
* Compare results against published SNS-TFQKD literature.

---

# Current Status Summary

### Completed

* Literature review of SNS-TFQKD.
* Python-based SNS-TFQKD simulation framework.
* Initial protocol performance evaluation.
* Analysis of provided OptiSystem architecture.
* Study of OptiSystem and MATLAB integration.

### Ongoing

* SNS-TFQKD protocol verification.
* OptiSystem model assessment.
* Send/Not-Send implementation development.

### Planned

* Decoy-state implementation.
* Detection statistics extraction.
* Secret key rate estimation.
* FSO channel integration.
* Comprehensive protocol validation.
