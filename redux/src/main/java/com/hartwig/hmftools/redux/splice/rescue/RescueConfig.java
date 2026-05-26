package com.hartwig.hmftools.redux.splice.rescue;

// Config for JunctionRescueResolver. Defaults mirror STAR's annotated-junction policy:
// alignSJoverhangMin=3, alignIntronMin=21, alignIntronMax=1_000_000.
public class RescueConfig
{
    public static final int DEFAULT_MIN_INTRON_LENGTH = 21;
    public static final int DEFAULT_MAX_INTRON_LENGTH = 1_000_000;
    public static final int DEFAULT_MIN_ANCHOR_OVERHANG = 3;
    public static final int DEFAULT_MAX_CHAIN_DEPTH = 4;
    // Tolerance (in read bases) for overlap between primary's matched region and supp's matched
    // region. BWA's primary and supplementary often disagree by 1-5 bases on where the junction
    // sits because they use different seed/extend choices. Allowing a small overlap and snapping
    // to an annotated junction within the disputed window catches the real splice without merging
    // garbage. STAR-style annotated-junction snap is preferred; ties resolve via "trust primary".
    public static final int DEFAULT_SOFTCLIP_TOLERANCE = 5;

    public final boolean Enabled;
    public final int MinIntronLength;
    public final int MaxIntronLength;
    public final int MinAnchorOverhang;
    public final int MaxChainDepth;
    public final boolean AnnotatedOnly;
    public final int SoftclipTolerance;

    public RescueConfig(
            final boolean enabled, final int minIntronLength, final int maxIntronLength,
            final int minAnchorOverhang, final int maxChainDepth, final boolean annotatedOnly,
            final int softclipTolerance)
    {
        Enabled = enabled;
        MinIntronLength = minIntronLength;
        MaxIntronLength = maxIntronLength;
        MinAnchorOverhang = minAnchorOverhang;
        MaxChainDepth = maxChainDepth;
        AnnotatedOnly = annotatedOnly;
        SoftclipTolerance = softclipTolerance;
    }

    public static RescueConfig defaults()
    {
        // tolerance=0 preserves strict-complementary semantics for tests that pre-date the snap
        // logic; production paths use enabledDefaults().
        return new RescueConfig(
                false,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MIN_ANCHOR_OVERHANG, DEFAULT_MAX_CHAIN_DEPTH,
                true, 0);
    }

    public static RescueConfig enabledDefaults()
    {
        // AnnotatedOnly=false in production: the merge constructs a better alignment from existing
        // BWA records (primary + supp), so annotation absence shouldn't block it. When annotation
        // does match a candidate snap point we still prefer that one (lets us recover STAR's exact
        // junction when annotated); otherwise trust-primary as fallback.
        return new RescueConfig(
                true,
                DEFAULT_MIN_INTRON_LENGTH, DEFAULT_MAX_INTRON_LENGTH,
                DEFAULT_MIN_ANCHOR_OVERHANG, DEFAULT_MAX_CHAIN_DEPTH,
                false, DEFAULT_SOFTCLIP_TOLERANCE);
    }
}
