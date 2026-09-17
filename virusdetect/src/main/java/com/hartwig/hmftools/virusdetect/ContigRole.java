package com.hartwig.hmftools.virusdetect;

// A candidate viral contig's standing within its oncology group after representative contig selection.
// "Abundant" means its vote share is near the oncology group's top; other candidates are low-abundance.
public enum ContigRole
{
    // The chosen representative: an abundant contig no abundant peer challenges, with the most read votes.
    REPRESENTATIVE,
    // An abundant contig indistinguishable from the representative (neither challenges the other); lost only the vote tie-break.
    REPRESENTATIVE_TWIN,
    // An abundant contig that an abundant peer challenges, so not the lead; it may still challenge others.
    SECONDARY,
    // An abundant contig unchallenged by peers but not crowned, its oncology group being unresolved.
    CONTESTED,
    // A low-abundance contig that challenges no abundant contig; a trace presence.
    MINOR,
    // A low-abundance contig that challenges an abundant contig, above its abundance; a possible hidden strain.
    MINOR_CHALLENGER
}
