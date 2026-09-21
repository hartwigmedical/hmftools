package com.hartwig.hmftools.viridian.selection;

import java.util.Set;

import com.hartwig.hmftools.viridian.detection.contig_support.ContigSupport;
import com.hartwig.hmftools.viridian.reference.ViralContig;

// Representative selection's verdict on one candidate contig, with the standing among its peers it was based on.
public record RepresentativeContigCandidate(
        ContigSupport support,
        // Position among the oncology group's candidates by read votes. 1 = best.
        int votesRank,
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
