package com.hartwig.hmftools.redux.consensus;

import htsjdk.samtools.SAMRecord;

public class ConsensusReadInfo
{
    public final SAMRecord ConsensusRead;

    // the actual read on which the consensus is based for all non-consensus properties
    public final SAMRecord TemplateRead;

    public final ConsensusOutcome Outcome;

    public ConsensusReadInfo(final SAMRecord consensusRead, final SAMRecord templateRead, final ConsensusOutcome outcome)
    {
        ConsensusRead = consensusRead;
        TemplateRead = templateRead;
        Outcome = outcome;
    }
}
