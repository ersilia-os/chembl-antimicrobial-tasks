"""Single source of truth for the stage-4 figures' semantic colors.

Every stage-4 figure (steps 21-27) pulls its colors from here, so the palette can never
drift across figures. Colors anchor to stylia's ``ArticleColors`` — the non-branded NPG
(Nature Publishing Group) publication palette used when the figures render under
``stylia.set_style("article")`` (set in ``plot_base`` / ``plot_utils``).

Design pattern mirrors ``zairachem/report/colors.py``: a semantic key -> stylia colour
name map is THE authority for what each colour means. Note that ``article`` has a
*different* set of named colours than ``ersilia`` (crimson/cobalt/tangerine/... rather
than plum/blue/mint/...), which is why the stage-4 scripts needed remapping to move off
the ersilia palette.

Conventions (used consistently across every stage-4 figure):
  - actives = crimson, inactives = silver.
  - Dose Response (DR) = cobalt, Single Point (SP) = tangerine.

Depends only on ``stylia`` so it never creates an import cycle with ``plot_base``.
"""

from stylia.colors import ArticleColors, CategoricalPalette, CyclicColormap, FadingColormap

_AC = ArticleColors()

# Semantic key -> stylia ArticleColors name. THE source of truth for colour meaning.
_SEMANTIC = {
    # activity outcome (binarisation figures)
    "active": "crimson",
    "inactive": "silver",
    # measurement format (used everywhere DR / SP appear)
    "DR": "cobalt",
    "SP": "tangerine",
    # eos embedding category split (step 22)
    "other": "silver",
    "dr_only": "cobalt",
    "sp_only": "tangerine",
    "shared": "orchid",
    # curation attrition (step 21 waterfall)
    "survivor": "cobalt",
    "lost": "amber",
    # class-composition stack (step 22)
    "all-inactive": "silver",
    "mixed": "amber",
    "all-active": "crimson",
    # per-panel decorative accents (step 21 overview; one colour per metric panel)
    "datasets": "cobalt",
    "chem_space": "crimson",
    "molecules": "periwinkle",
    "strains": "turquoise",
    "sizes": "amber",
    "reuse": "orchid",
    "docs": "tangerine",
    # target-type regions (step 27 curation story)
    "wholecell": "turquoise",
    "protein": "amber",
    "wholecell_protein": "orchid",
    "other_region": "silver",
    # cross-pathogen curation accents (step 27)
    "total": "cobalt",          # "% of total chemical space"
    "related": "crimson",       # ChEMBL-coverage highlight slice
    # generic
    "neutral": "silver",
    "ink": "black",             # structural box/whisker lines
}


def color(key):
    """RGB tuple for a semantic key (for matplotlib)."""
    return getattr(_AC, _SEMANTIC[key])


def named(name):
    """Direct access to a stylia article colour by its own name — for the fixed
    multi-category colour lists below where a per-key semantic name adds no clarity."""
    return getattr(_AC, name)


def lighten(key, amount):
    """Lightened RGB for a semantic key (e.g. a faint 'lost' cap over a solid bar)."""
    return _AC.get(_SEMANTIC[key], lighten=amount)


def lighten_named(name, amount):
    """Lightened RGB for a stylia article colour by its own name (faded-vs-full pairs, e.g.
    a pre/post 'before' box drawn in a lighter shade of the 'after' colour)."""
    return _AC.get(name, lighten=amount)


def category_palette(n):
    """``n`` perceptually distinct colours for a categorical series (NPG family)."""
    return CategoricalPalette("npg").get(max(n, 1))


def cyclic(n):
    """``n`` colours from the cyclic NPG colormap (used for the top-pool clusters and
    the ordered-bin stacks: size / hit-count / compatibility-group bins)."""
    cm = CyclicColormap("npg")
    cm.fit(list(range(n + 1)))
    return cm.transform(list(range(n)))


# --- fixed ordered colour lists for step-23 multi-category stacks ------------
# Kept explicit (rather than a palette slice) so each category keeps a stable,
# legible colour across every cross-pathogen figure.

def targettype_colors():
    """TARGET_TYPES order: ORGANISM (whole-cell), SINGLE PROTEIN, DISCARDED."""
    return [color("wholecell"), color("protein"), color("other_region")]


def region_colors():
    """TT_MOL_REGIONS order: Whole-cell, Whole-cell+Protein, Protein target, Other."""
    return [color("wholecell"), color("wholecell_protein"),
            color("protein"), color("other_region")]


def composition_colors():
    """COMPOSITION order: all-inactive, mixed, all-active."""
    return [color("all-inactive"), color("mixed"), color("all-active")]


def reuse_colors():
    """Reuse bins: 1 dataset, 2 datasets, >=3 datasets (low -> high redundancy)."""
    return [named("silver"), named("cobalt"), named("crimson")]


def discard_step_colors():
    """STEP_ORDER order (9): retained, not-whole-cell, qualitative-only, uncategorized,
    Co-ADD, PubChem, <=5 molecules, multi-unit, no-direction. Retained kept calm (cobalt);
    reasons distinct. Must stay length- and order-aligned with plot_utils.STEP_ORDER."""
    return [named("cobalt"), named("amber"), named("tangerine"), named("turquoise"),
            named("orchid"), named("periwinkle"), named("fuchsia"), named("silver"),
            named("crimson")]


def fading(key="total"):
    """A sequential (fading) colormap in an article hue, for magnitude heatmaps."""
    hue = {"total": "cobalt", "related": "crimson", "active": "crimson"}.get(key, "cobalt")
    return FadingColormap(hue)
