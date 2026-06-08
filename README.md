Singing Voice Dataset
This repository contains a singing voice dataset designed for the acoustic analysis of vocal developmental stages and the effects of vocal pedagogical intervention. The dataset consists of short audio fragments of the word "tanto" from the Italian art song <Caro mio ben> composed by T. Giordani.
1. Recording and Acoustic Analysis Conditions
•	Equipment & Environment: Olympus LS-P2 IC recorder / Vocal lesson room (Microphone-to-subject distance: approx. 2 m)
•	Audio Format: 44.1 kHz, 16-bit, Mono
•	Acoustic Analysis Parameters: Analysis window length (frame size) of 2,048 samples
2. Subjects and Dataset Structure All subjects are soprano singers. For each student subject, audio data was recorded at two distinct stages: Z1 (at admission) and Z2 (after two years of vocal instruction). Subject Groups Group ID N Description Education E01–E07 7 Students in the Music Education Program Music College V01–V13 13 Students majoring in Vocal Performance Reference Model M2025 1 Professional soprano (normative pedagogical model for observing developmental progress)
3. Naming Convention
[Subject ID][Stage]_[Detailed Information].wav
•	Stage Mapping: before = Z1 (At admission) / after = Z2 (After 2 years of instruction)
•	Detailed Information: Includes the recording date, song fragment identifier (tanto), and channel type (mono).
Examples
•	E01before_tanto-tanto_20141001_02_mono.wav
o	Audio data for Subject E01 (Education) at admission (Z1). Recorded on October 1, 2014.
•	V05after_20240502_tanto_tanto_02_mono.wav
o	Audio data for Subject V05 (Music College) after two years of instruction (Z2). Recorded on May 2, 2024.
•	M2025.wav
o	Pedagogical reference model recorded by a professional singer.
