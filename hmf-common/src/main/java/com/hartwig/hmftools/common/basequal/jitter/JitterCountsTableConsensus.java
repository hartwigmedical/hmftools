package com.hartwig.hmftools.common.basequal.jitter;

import com.hartwig.hmftools.common.qual.BqrReadType;

public class JitterCountsTableConsensus extends JitterCountsTable
{
    public final BqrReadType ConsensusType;

    public JitterCountsTableConsensus(
            final String repeatUnit, final double maxSingleAltSiteContributionPerc, final BqrReadType consensusType)
    {
        super(repeatUnit, maxSingleAltSiteContributionPerc);
        ConsensusType = consensusType;
    }
}