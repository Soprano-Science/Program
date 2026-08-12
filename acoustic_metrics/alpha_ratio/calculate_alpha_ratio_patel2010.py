from __future__ import annotations

import platform
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import scipy
from scipy.signal import welch
import soundfile as sf


WAV_FILE = Path(
    "kawano-after_2022_20190408tanto_tanto_02_mono(3).wav"
)

# The uploaded file contains the word "tanto". Patel et al. (2010)
# analyzed vowel /a/, so this script uses the same stable /a/ interval
# used in the preceding Q-value analysis.
START_SEC = 0.430
END_SEC = 0.930

# Reproducible LTAS implementation. Patel et al. do not report the
# internal FFT/window settings of the Cofi software.
N_FFT = 2048
N_OVERLAP = N_FFT // 2

LOW_MIN_HZ = 50.0
LOW_MAX_HZ = 1000.0
HIGH_MIN_HZ = 1000.0
HIGH_MAX_HZ = 5000.0


def package_version(package_name: str) -> str:
    try:
        return version(package_name)
    except PackageNotFoundError:
        return "not installed"


def calculate_alpha_ratio(
    audio: np.ndarray,
    sample_rate: int,
) -> dict[str, float]:
    
"""Calculate the alpha-ratio measure from an LTAS estimate using the
frequency bands described by Patel et al. (2010).

Alpha ratio was defined in this study as:

    alpha_linear = E(50–1000 Hz) / E(1–5 kHz)
    alpha_dB = 10 * log10(alpha_linear)

Thus, the orientation used in this study is low-frequency energy
divided by high-frequency energy. Accordingly, lower alpha-ratio
values indicate relatively greater high-frequency energy in the
1–5-kHz band.

Because the direction of the alpha-ratio convention is not always
stated consistently in the literature, the numerator and denominator
are specified explicitly here. The measure is interpreted
descriptively as a measure of spectral energy balance and not as a
direct physiological measure.
"""

    frequencies, psd = welch(
        audio,
        fs=sample_rate,
        window="hann",
        nperseg=N_FFT,
        noverlap=N_OVERLAP,
        nfft=N_FFT,
        detrend=False,
        scaling="density",
        return_onesided=True,
    )

    frequency_step = frequencies[1] - frequencies[0]

    low_mask = (
        (frequencies >= LOW_MIN_HZ)
        & (frequencies < LOW_MAX_HZ)
    )
    high_mask = (
        (frequencies >= HIGH_MIN_HZ)
        & (frequencies <= HIGH_MAX_HZ)
    )

    # Integrate the PSD in each band. Multiplication by df makes these
    # band-energy estimates; df cancels in the ratio but is retained for
    # clarity.
    low_energy = float(np.sum(psd[low_mask]) * frequency_step)
    high_energy = float(np.sum(psd[high_mask]) * frequency_step)

    if low_energy <= 0.0 or high_energy <= 0.0:
        raise ValueError("Band energy must be positive.")

    alpha_linear = high_energy / low_energy
    alpha_db = 10.0 * np.log10(alpha_linear)

    return {
        "frequency_step_hz": float(frequency_step),
        "low_energy": low_energy,
        "high_energy": high_energy,
        "alpha_linear": float(alpha_linear),
        "alpha_db": float(alpha_db),
        "inverse_linear": float(1.0 / alpha_linear),
        "inverse_db": float(-alpha_db),
    }


def main() -> None:
    audio, sample_rate = sf.read(WAV_FILE, always_2d=False)

    if audio.ndim != 1:
        raise ValueError("A mono WAV file is required.")

    start_sample = round(START_SEC * sample_rate)
    end_sample = round(END_SEC * sample_rate)

    if start_sample < 0 or end_sample > len(audio):
        raise ValueError("The selected interval is outside the audio file.")

    segment = np.asarray(audio[start_sample:end_sample], dtype=np.float64)
    result = calculate_alpha_ratio(segment, sample_rate)

    print("[Environment]")
    print(f"Python:    {platform.python_version()}")
    print(f"NumPy:     {np.__version__}")
    print(f"SciPy:     {scipy.__version__}")
    print(f"SoundFile: {sf.__version__}")
    print(f"Librosa:   {package_version('librosa')} (not used in calculation)")
    print()

    print("[Analysis interval]")
    print(f"Sample rate: {sample_rate} Hz")
    print(f"Start:       {start_sample / sample_rate:.6f} s")
    print(f"End:         {end_sample / sample_rate:.6f} s")
    print(f"Duration:    {len(segment) / sample_rate:.6f} s")
    print(f"FFT length:  {N_FFT}")
    print(f"Overlap:     {N_OVERLAP}")
    print(f"Frequency resolution: {result['frequency_step_hz']:.12f} Hz")
    print()

    print("[Band energies]")
    print(f"E(50-1000 Hz): {result['low_energy']:.15f}")
    print(f"E(1-5 kHz):    {result['high_energy']:.15f}")
    print()

    print("[Alpha ratio]")
    print(f"Linear E_high/E_low: {result['alpha_linear']:.15f}")
    print(f"Alpha ratio:         {result['alpha_db']:.12f} dB")
    print()
    print("For the inverse low/high convention:")
    print(f"Linear E_low/E_high: {result['inverse_linear']:.15f}")
    print(f"Inverse dB value:    {result['inverse_db']:.12f} dB")


if __name__ == "__main__":
    main()
