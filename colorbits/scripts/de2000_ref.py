"""Independent CIEDE2000 reference, transcribed line-by-line from
Sharma, Wu & Dalal (2005), Color Research and Application 30(1), Table I.
Vectorized for (..., 3) Lab arrays."""
import numpy as np


def _hp_deg(bp, ap):
    # arctan2 gives (-180, 180]; mod 360 maps negative angles into [0, 360)
    return np.degrees(np.arctan2(bp, ap)) % 360.0


def de2000_ref(lab1, lab2, kL=1.0, kC=1.0, kH=1.0):
    lab1 = np.asarray(lab1, dtype=np.float64)
    lab2 = np.asarray(lab2, dtype=np.float64)
    L1, a1, b1 = lab1[..., 0], lab1[..., 1], lab1[..., 2]
    L2, a2, b2 = lab2[..., 0], lab2[..., 1], lab2[..., 2]

    C1 = np.hypot(a1, b1)
    C2 = np.hypot(a2, b2)
    Cbar = 0.5 * (C1 + C2)
    G = 0.5 * (1 - np.sqrt(Cbar**7 / (Cbar**7 + 25.0**7)))
    ap1 = (1 + G) * a1
    ap2 = (1 + G) * a2
    Cp1 = np.hypot(ap1, b1)
    Cp2 = np.hypot(ap2, b2)
    hp1 = _hp_deg(b1, ap1)
    hp2 = _hp_deg(b2, ap2)

    dLp = L2 - L1
    dCp = Cp2 - Cp1
    dhp = np.where(
        Cp1 * Cp2 == 0,
        0.0,
        np.where(
            np.abs(hp2 - hp1) <= 180,
            hp2 - hp1,
            np.where(hp2 - hp1 > 180, hp2 - hp1 - 360, hp2 - hp1 + 360)
        ),
    )
    dHp = 2 * np.sqrt(Cp1 * Cp2) * np.sin(np.radians(dhp / 2.0))

    Lbp = (L1 + L2) / 2.0
    Cbp = (Cp1 + Cp2) / 2.0
    hbp = np.where(
        Cp1 * Cp2 == 0,
        hp1 + hp2,
        np.where(
            np.abs(hp1 - hp2) <= 180,
            (hp1 + hp2) / 2.0,
            np.where(hp1 + hp2 < 360, (hp1 + hp2 + 360) / 2.0,
                     (hp1 + hp2 - 360) / 2.0),
        ),
    )
    T = (1
         - 0.17 * np.cos(np.radians(hbp - 30))
         + 0.24 * np.cos(np.radians(2 * hbp))
         + 0.32 * np.cos(np.radians(3 * hbp + 6))
         - 0.20 * np.cos(np.radians(4 * hbp - 63)))
    dtheta = 30 * np.exp(-(((hbp - 275) / 25.0) ** 2))
    RC = 2 * np.sqrt(Cbp**7 / (Cbp**7 + 25.0**7))
    SL = 1 + (0.015 * (Lbp - 50) ** 2) / np.sqrt(20 + (Lbp - 50) ** 2)
    SC = 1 + 0.045 * Cbp
    SH = 1 + 0.015 * Cbp * T
    RT = -np.sin(np.radians(2 * dtheta)) * RC

    return np.sqrt(
        (dLp / (kL * SL)) ** 2
        + (dCp / (kC * SC)) ** 2
        + (dHp / (kH * SH)) ** 2
        + RT * (dCp / (kC * SC)) * (dHp / (kH * SH))
    )
