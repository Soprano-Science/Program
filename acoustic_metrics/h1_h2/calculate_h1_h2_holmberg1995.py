import numpy as np
import soundfile as sf
import parselmouth
from scipy import signal

WAV_FILE = (
    "kawano-after_2022_20190408tanto_"
    "tanto_02_mono(3).wav"
)

CENTER_SEC = 0.680
TARGET_SR = 10000
FRAME_LENGTH = 512       # 51.2 ms
LOWPASS_HZ = 4500.0
N_FFT = 131072


def correction_db(harmonic_hz, f1_hz):
    """Return the correction term for the influence of F1 on harmonic amplitude."""
    denominator = abs(
        1.0 - (harmonic_hz / f1_hz) ** 2
    )
    return 20.0 * np.log10(1.0 / denominator)


# Load the audio file
audio, original_sr = sf.read(
    WAV_FILE,
    always_2d=False,
)
audio = np.asarray(audio, dtype=np.float64)

# Apply a 4.5 kHz low-pass filter
sos = signal.butter(
    8,
    LOWPASS_HZ,
    btype="lowpass",
    fs=original_sr,
    output="sos",
)
filtered = signal.sosfiltfilt(sos, audio)

# Resample from 44.1 kHz to 10 kHz
audio_10k = signal.resample_poly(
    filtered,
    up=100,
    down=441,
)

# Extract a 51.2 ms frame centered on the vowel midpoint
center_sample = round(CENTER_SEC * TARGET_SR)
start = center_sample - FRAME_LENGTH // 2
frame = audio_10k[start:start + FRAME_LENGTH]

frame = frame - np.mean(frame)
frame = frame * np.hamming(FRAME_LENGTH)

# Calculate the amplitude spectrum
spectrum = np.fft.rfft(frame, n=N_FFT)
frequencies = np.fft.rfftfreq(
    N_FFT,
    d=1.0 / TARGET_SR,
)
spectrum_db = 20.0 * np.log10(
    np.maximum(np.abs(spectrum), 1e-15)
)

# Estimate F0 using Praat
sound = parselmouth.Sound(WAV_FILE)

pitch = sound.to_pitch_ac(
    time_step=0.001,
    pitch_floor=400.0,
    pitch_ceiling=800.0,
)
f0 = pitch.get_value_at_time(CENTER_SEC)

# Detect spectral peaks
peak_indices, _ = signal.find_peaks(spectrum_db)

# H1: select the largest peak near F0
h1_candidates = peak_indices[
    (frequencies[peak_indices] >= 0.8 * f0)
    & (frequencies[peak_indices] <= 1.2 * f0)
]
h1_index = h1_candidates[
    np.argmax(spectrum_db[h1_candidates])
]

# H2: select the largest peak near 2 x F0
h2_candidates = peak_indices[
    (frequencies[peak_indices] >= 1.6 * f0)
    & (frequencies[peak_indices] <= 2.4 * f0)
]
h2_index = h2_candidates[
    np.argmax(spectrum_db[h2_candidates])
]

h1_hz = frequencies[h1_index]
h2_hz = frequencies[h2_index]

h1_db = spectrum_db[h1_index]
h2_db = spectrum_db[h2_index]

raw_h1_h2 = h1_db - h2_db

# Estimate F1 using Praat Formant (Burg)
formant = sound.to_formant_burg(
    time_step=0.001,
    max_number_of_formants=5,
    maximum_formant=5000.0,
    window_length=0.0512,
    pre_emphasis_from=50.0,
)
f1_hz = formant.get_value_at_time(
    1,
    CENTER_SEC,
)

# Apply the F1 correction
h1_corrected = (
    h1_db - correction_db(h1_hz, f1_hz)
)
h2_corrected = (
    h2_db - correction_db(h2_hz, f1_hz)
)

corrected_h1_h2 = (
    h1_corrected - h2_corrected
)

print(f"H1: {h1_hz:.6f} Hz")
print(f"H2: {h2_hz:.6f} Hz")
print(f"F1: {f1_hz:.6f} Hz")
print(f"Uncorrected H1-H2: {raw_h1_h2:.6f} dB")
print(
    f"F1-corrected H1*-H2*: "
    f"{corrected_h1_h2:.6f} dB"
)
