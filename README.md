<h1>Singing Voice Dataset</h1>

<p>This repository contains a singing voice dataset designed for the acoustic analysis of vocal developmental stages and the effects of vocal pedagogical intervention. The dataset consists of short audio fragments of the word &quot;tanto&quot; from the Italian art song &lt;<em>Caro mio ben</em>&gt; composed by T. Giordani.</p>

<hr>

<h2>1. Recording and Acoustic Analysis Conditions</h2>

<ul>
  <li><strong>Equipment &amp; Environment</strong>: Olympus LS-P2 IC recorder / Vocal lesson room (Microphone-to-subject distance: approx. 2 m)</li>
  <li><strong>Audio Format</strong>: 44.1 kHz, 16-bit, Mono</li>
  <li><strong>Acoustic Analysis Parameters</strong>: Analysis window length (frame size) of 2,048 samples</li>
</ul>

<hr>
2. Subjects and Dataset Structure

All subjects are **soprano** singers. For each student subject, audio data was recorded at two distinct stages: **Z1 (at admission)** and **Z2 (after two years of vocal instruction)**.

<hr>
<p><strong>Subject Groups</strong></p>

<table border="0" cellpadding="6" cellspacing="0" style="border-collapse: collapse; width: 100%;">
  <thead>
    <tr style="border-top: 2px solid #000; border-bottom: 1px solid #000; font-weight: bold; text-align: left;">
      <th style="padding: 8px;">Group</th>
      <th style="padding: 8px;">ID</th>
      <th style="padding: 8px;">N</th>
      <th style="padding: 8px;">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 8px;"><strong>Education</strong></td>
      <td style="padding: 8px;">E01–E07</td>
      <td style="padding: 8px;">7</td>
      <td style="padding: 8px;">Students in the Music Education Program</td>
    </tr>
    <tr>
      <td style="padding: 8px;"><strong>Music College</strong></td>
      <td style="padding: 8px;">V01–V13</td>
      <td style="padding: 8px;">13</td>
      <td style="padding: 8px;">Students majoring in Vocal Performance</td>
    </tr>
    <tr style="border-bottom: 2px solid #000;">
      <td style="padding: 8px;"><strong>Reference Model</strong></td>
      <td style="padding: 8px;">M2025</td>
      <td style="padding: 8px;">1</td>
      <td style="padding: 8px;">Professional soprano (normative pedagogical model for observing developmental progress)</td>
    </tr>
  </tbody>
</table>
<hr>
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
