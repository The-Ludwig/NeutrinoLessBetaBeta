import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# This did not work, maybe some https error?
# URL = "http://www-nds.iaea.org/amdc/ame2020/nubase_3.mas20.txt"
U_TO_KEV = 931494.10242
URL = "./data.txt"
NAMES = (
    "A",
    "Z",
    "i",
    "AEl",
    "s",
    "ME",
    "uME",
    "Exc",
    "dE",
    "Orig",
    "Isom.Unc",
    "Isom.Inv",
    "T",
    "unit T",
    "dT",
    "Jpi",
    "Ensdf",
    "Discovery",
    "BR",
)
COL_STARTS = (
    0,
    4,
    7,
    11,
    16,
    18,
    31,
    42,
    54,
    65,
    67,
    68,
    69,
    78,
    81,
    88,
    102,
    114,
    119,
)
COL_ENDS = (2, 6, 8, 15, 16, 30, 41, 53, 64, 66, 67, 68, 77, 79, 87, 101, 103, 117, 208)


def get_table():
    data = ascii.read(
        URL,
        format="fixed_width_no_header",
        col_starts=COL_STARTS,
        col_ends=COL_ENDS,
        names=NAMES,
    )

    data = data[~data["ME"].mask]

    data.add_column(
        [ael.replace(f"{a}", "") for ael, a in zip(data["AEl"], data["A"])],
        name="Symbol",
        index=2,
    )
    data.add_column(
        [
            np.ma.core.MaskedConstant
            if isinstance(me, np.ma.core.MaskedConstant)
            else float(me.replace("#", ""))
            for me in data["ME"]
        ],
        name="M-A",
        index=2,
    )
    data.add_column(
        [
            me + a.view(np.ma.MaskedArray) * U_TO_KEV
            for me, a in zip(data["M-A"], data["A"])
        ],
        name="M",
        index=3,
    )
    return data


def make_plot(
    filename,
    width=4,
    yval="M-A",
    A=76,
    beta=False,
    betabeta=False,
    element=True,
    max_m=800e6,
    betabeta_text=False,
    beta_arrow=True,
    title=False,
    scale=1e-3,
):
    golden_ratio = (5 ** (0.5) - 1) / 2
    print(f"size = {width}, {width*golden_ratio}")
    plt.gcf().set_size_inches(width, width * golden_ratio)

    data = get_table()
    A76 = data[
        np.all(np.array([data["A"] == A, data["i"] == 0, data[yval] < max_m]), axis=0)
    ]

    plt.plot(A76["Z"], A76[yval] * scale, "kx")
    if element:
        for el in A76:
            plt.text(el["Z"] + 0.1, (el[yval]) * scale, el["Symbol"])

    m_e = 511

    mask_sorted = np.argsort(A76["Z"])
    for par, daug, daugdaug in zip(
        A76[mask_sorted[:-2]], A76[mask_sorted[1:-1]], A76[mask_sorted[2:]]
    ):
        dm = daug[yval] - par[yval]
        ddm = daugdaug[yval] - par[yval]

        if -dm > m_e:
            if beta:
                plt.text(par["Z"] + 0.6, (par[yval] + 0.5 * dm) * scale, r"$\beta^-$")
            if beta_arrow:
                plt.annotate(
                    "",
                    xytext=(par["Z"], par[yval] * scale),
                    xy=(daug["Z"], daug[yval] * scale),
                    arrowprops=dict(arrowstyle="->", color="blue"),
                )
        elif -ddm > 2 * m_e and betabeta:
            if betabeta_text:
                plt.text(
                    par["Z"] + 0.6,
                    (par[yval] - 0.5 * ddm - 3000) * scale,
                    r"$\beta^- \beta^-$",
                    va="top",
                )
            plt.annotate(
                "",
                xytext=(par["Z"], par[yval] * scale),
                xy=(daugdaug["Z"], daugdaug[yval] * scale),
                arrowprops=dict(arrowstyle="->", color="blue"),
            )

    if title:
        plt.title(f"$A={A}$")
    plt.xlabel("Protons $Z$")
    plt.ylabel(r"Binding Energy $-E / \si{\mega\electronvolt}$")

    plt.xlim(min(A76["Z"]) - 1, max(A76["Z"]) + 1)
    plt.ylim((min(A76[yval]) - 5000) * scale, (max(A76[yval]) + 5000) * scale)

    plt.tight_layout()
    plt.savefig(filename, transparent=True)
    plt.cla()


if __name__ == "__main__":
    make_plot("build/const_a_all.pgf")
    make_plot("build/const_a_zoom.pgf", max_m=-60e3, beta=True)
    make_plot(
        "build/const_a_betabeta.pgf",
        max_m=-60e3,
        betabeta=True,
        betabeta_text=True,
        beta_arrow=False,
    )
    make_plot(
        "build/const_a_XE.pgf",
        A=136,
        betabeta=True,
        betabeta_text=False,
        beta_arrow=False,
    )
