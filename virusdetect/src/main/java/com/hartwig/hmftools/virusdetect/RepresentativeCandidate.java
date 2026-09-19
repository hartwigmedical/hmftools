package com.hartwig.hmftools.virusdetect;

import java.util.Set;

// Representative selection's verdict on one candidate contig, with the standing among its peers it was based on.
public record RepresentativeCandidate(
        ContigSupport support,
        // Read votes near the group's highest, so abundant enough to contest the lead.
        boolean comparable,
        // Peers this contig fits a decisive share of the group's reads better than, and those which do so to it.
        Set<ViralContig> challenges,
        Set<ViralContig> challengedBy,
        ContigRole role
)
{
    public ViralContig contig()
    {
        return support.contig();
    }
}
