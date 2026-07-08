"""Publication-figure framework for the stage-4 plotting scripts.

Adopts the organisation of ``zairachem/report`` in a form that fits this repo's
one-pathogen-per-invocation scripts:

  - Applies the publication house style (``print`` format + ``article`` style) at
    import, so any script that imports this module renders paper-ready figures.
    ``plot_utils`` applies the same style at its own import, so the scripts that use
    ``plot_utils`` but not ``plot_base`` (steps 23-27) stay consistent.
  - ``BasePlot``: each concrete panel is a subclass that draws into a *provided* stylia
    axis (never creating its own figure — per the stylia house rules) and records the
    tidy data it plotted in ``self.data`` (or ``None``). An ``is_available`` guard lets a
    panel render a labelled blank instead of crashing on empty input.
  - ``render_figure``: builds a stylia grid, hands each panel its axis, optionally places
    a shared legend *outside* the data area, and saves. Multi-panel grids use stylia's
    default width/height (``save_figure`` uses ``bbox_inches="tight"``, so rotated tick
    labels are never clipped).
  - Legend helpers (``figure_legend`` / ``axis_legend``): always frameless and positioned
    outside the plotting area — never ``loc="best"``.
  - ``EmbeddingFigure``: the 2x2 four-projection eos1klk scatter, article-coloured, with a
    single cleaned bottom legend.
"""

import os

import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

import stylia


def _compact(x, _pos=None):
    """Human-readable large-number tick label: 200000 -> '200k', 1.5e6 -> '1.5M'."""
    a = abs(x)
    if a >= 1e6:
        return f"{x / 1e6:g}M"
    if a >= 1e3:
        return f"{x / 1e3:g}k"
    return f"{x:g}"


def compact_yaxis(ax):
    """Format an axis' y ticks compactly (k / M) so large counts stay short."""
    ax.yaxis.set_major_formatter(FuncFormatter(_compact))

# Format: print | Style: article - publication figures. Applied at import so any importer
# renders paper-ready; plot_utils sets the same style for the scripts that skip plot_base.
stylia.set_format("print")
stylia.set_style("article")


# --- base panel -------------------------------------------------------------

class BasePlot:
    """A single figure panel. Subclasses draw into ``ax`` in their ``__init__`` and set
    ``self.data`` (a tidy DataFrame of exactly what was plotted, or None)."""

    name = "base"

    def __init__(self, ax, abc=None):
        self.ax = ax
        self.abc = abc
        self.data = None
        self.is_available = True

    def _label(self, xlabel="", ylabel="", title=None):
        stylia.label(self.ax, xlabel=xlabel, ylabel=ylabel, title=title, abc=self.abc)

    def _empty(self, title=None, xlabel="", ylabel=""):
        """Render a labelled blank panel (no data available)."""
        self.is_available = False
        self._label(xlabel=xlabel, ylabel=ylabel, title=title)


# --- legend helpers (always outside the data, frameless) --------------------

def _handles(entries):
    return [Patch(color=c, label=l) for l, c in entries]


def figure_legend(fig, entries, loc="lower center", anchor=(0.5, -0.01), ncol=None, title=None):
    """Shared, figure-level legend placed outside the panels (default: below them)."""
    fig.legend(handles=_handles(entries), loc=loc, ncol=ncol or len(entries),
               bbox_to_anchor=anchor, frameon=False, title=title)


def axis_legend(ax, entries, loc="upper left", anchor=(1.02, 1.0), ncol=1, title=None, frame=False):
    """Legend anchored just outside a single axis (default: to its right).

    ``frame=True`` draws a semi-transparent white background (no border) so an in-axis
    legend stays readable over busy data; default stays frameless."""
    leg = ax.legend(handles=_handles(entries), loc=loc, bbox_to_anchor=anchor,
                    ncol=ncol, frameon=frame, title=title)
    if frame and leg is not None:
        leg.get_frame().set_facecolor("white")
        leg.get_frame().set_alpha(0.7)
        leg.get_frame().set_edgecolor("none")
    return leg


def inset_legend(ax, loc="upper right", title=None, ncol=1, entries=None):
    """In-axis legend with a semi-transparent white background so it stays readable over
    the data. Uses the axis' own labelled artists by default; pass ``entries`` (list of
    ``(label, colour)``) for artists that carry no label (e.g. boxplots). Use when the
    legend should sit inside the panel (its meaning already implied by the panel title)."""
    if entries is not None:
        leg = ax.legend(handles=_handles(entries), loc=loc, frameon=True, title=title, ncol=ncol)
    else:
        leg = ax.legend(loc=loc, frameon=True, title=title, ncol=ncol)
    if leg is not None:
        leg.get_frame().set_facecolor("white")
        leg.get_frame().set_alpha(0.7)
        leg.get_frame().set_edgecolor("none")
    return leg


# --- composite figure driver ------------------------------------------------

def render_figure(out_dir, name, nrows, ncols, builders, legend=None,
                  width=None, height=None, width_ratios=None, height_ratios=None):
    """Create a stylia grid, draw each panel via its builder (``ax -> BasePlot``), place
    an optional shared legend, and save ``<out_dir>/<name>.png``.

    ``legend`` is an optional callable ``(fig, panels) -> None`` so a figure can place a
    shared legend (or several) exactly where it wants, outside the data. ``width_ratios`` /
    ``height_ratios`` set unequal column/row sizes (forwarded to stylia). Returns the list
    of instantiated panels (their ``.data`` is available for companion CSVs)."""
    kwargs = {}
    if width is not None:
        kwargs["width"] = width
    if height is not None:
        kwargs["height"] = height
    if width_ratios is not None:
        kwargs["width_ratios"] = width_ratios
    if height_ratios is not None:
        kwargs["height_ratios"] = height_ratios
    fig, axs = stylia.create_figure(nrows, ncols, **kwargs)
    panels = [b(axs.next()) for b in builders]
    if legend is not None:
        legend(fig, panels)
    path = os.path.join(out_dir, name + ".png")
    stylia.save_figure(path)
    print(f"  figure -> {path}")
    return panels


# --- eos1klk embedding figure (2x2 projections) -----------------------------

class EmbeddingFigure:
    """2x2 scatter of the four eos1klk projections (PCA / t-MAP / t-SNE / UMAP).

    ``groups`` is a list of ``(label, bool-mask, colour)`` drawn in order (first = drawn
    first / background). A single frameless legend is placed centred below the panels.
    """

    def __init__(self, cc, groups):
        self.cc = cc
        self.groups = groups

    def save(self, out_dir, name, legend_ncol=None, companion=False):
        from plot_utils import EMBEDDINGS  # data-column constant; module already loaded

        cc = self.cc
        fig, axs = stylia.create_figure(2, 2, width=0.8, height=0.8)
        for proj, xcol, ycol in EMBEDDINGS:
            ax = axs.next()
            for label, mask, color in self.groups:
                if mask.any():
                    ax.scatter(cc[xcol].values[mask], cc[ycol].values[mask], s=2,
                               alpha=0.4, color=color, edgecolors="none", rasterized=True)
            ax.set_xticks([])
            ax.set_yticks([])
            stylia.label(ax, xlabel="", ylabel="", title=proj)
        # Anchor the legend by its TOP just below the panels (loc="upper center") so it
        # grows downward and never overlaps the bottom row, however many rows it needs.
        # bbox_inches="tight" (stylia.save_figure) keeps it in frame.
        figure_legend(fig, [(l, c) for l, m, c in self.groups],
                      loc="upper center", anchor=(0.5, -0.03),
                      ncol=legend_ncol or len(self.groups))
        path = os.path.join(out_dir, name + ".png")
        stylia.save_figure(path)
        print(f"  figure -> {path}")
        if companion:
            self._write_companion(out_dir, name, EMBEDDINGS)

    def _write_companion(self, out_dir, name, embeddings):
        """Tidy CSV of the plotted points: inchikey, every projection x/y, and group.
        These coordinates/memberships live in no existing stage-4 CSV."""
        cc = self.cc
        group_of = pd.Series("", index=cc.index)
        for label, mask, color in self.groups:  # last matching group wins (draw order)
            group_of.loc[cc.index[mask]] = label
        cols = {"inchikey": cc["inchikey"].values} if "inchikey" in cc.columns else {}
        for proj, xcol, ycol in embeddings:
            cols[f"{proj}_x"] = cc[xcol].values
            cols[f"{proj}_y"] = cc[ycol].values
        cols["group"] = group_of.values
        csv_path = os.path.join(out_dir, name + ".csv")
        pd.DataFrame(cols).to_csv(csv_path, index=False)
        print(f"  data   -> {csv_path}")
