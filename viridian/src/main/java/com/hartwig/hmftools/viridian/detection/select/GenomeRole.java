package com.hartwig.hmftools.viridian.detection.select;

// A candidate virus genome's standing within its oncology group after representative selection.
// "Abundant" means its read vote share is near the highest for that oncology group.
public enum GenomeRole
{
    // The chosen representative: an abundant genome no abundant peer challenges, with the most read votes.
    REPRESENTATIVE,
    // An abundant genome indistinguishable from the representative (neither challenges the other). Lost only
    // the vote tie-break.
    REPRESENTATIVE_TWIN,
    // An abundant genome that an abundant peer challenges, so not the lead. It may still challenge others.
    SECONDARY,
    // An abundant genome unchallenged by peers but not crowned, its oncology group being unresolved.
    CONTESTED,
    // A low-abundance genome that challenges no abundant genome. A trace presence.
    MINOR,
    // A low-abundance genome that challenges an abundant one, above its abundance. A possible hidden strain.
    MINOR_CHALLENGER
}
