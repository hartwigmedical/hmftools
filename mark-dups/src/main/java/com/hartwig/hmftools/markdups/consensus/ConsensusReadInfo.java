package com.hartwig.hmftools.markdups.consensus;

import htsjdk.samtools.SAMRecord;

public class ConsensusReadInfo
{
    public final SAMRecord ConsensusRead;
    public final ConsensusOutcome Outcome;

    public ConsensusReadInfo(final SAMRecord consensusRead, final ConsensusOutcome outcome)
    {
        ConsensusRead = consensusRead;
        Outcome = outcome;
    }
}
