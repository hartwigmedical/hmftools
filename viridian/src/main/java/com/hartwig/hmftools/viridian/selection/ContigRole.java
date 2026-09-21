package com.hartwig.hmftools.viridian.selection;

// A candidate viral contig's standing within its oncology group after representative contig selection.
// "Abundant" means its read vote share is near the highest for that oncology group.
public enum ContigRole
{
    // The chosen representative: an abundant contig no abundant peer challenges, with the most read votes.
    REPRESENTATIVE,
    // An abundant contig indistinguishable from the representative (neither challenges the other). Lost only the vote tie-break.
    REPRESENTATIVE_TWIN,
    // An abundant contig that an abundant peer challenges, so not the lead. It may still challenge others.
    SECONDARY,
    // An abundant contig unchallenged by peers but not crowned, its oncology group being unresolved.
    CONTESTED,
    // A low-abundance contig that challenges no abundant contig. A trace presence.
    MINOR,
    // A low-abundance contig that challenges an abundant contig, above its abundance. A possible hidden strain.
    MINOR_CHALLENGER
}
