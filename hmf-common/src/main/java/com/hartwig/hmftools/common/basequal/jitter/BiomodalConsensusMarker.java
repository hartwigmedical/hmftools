package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.basequal.jitter.ConsensusType.HIGH_QUAL;
import static com.hartwig.hmftools.common.basequal.jitter.ConsensusType.NONE;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.LOW_QUAL_CUTOFF;

import htsjdk.samtools.SAMRecord;

public class BiomodalConsensusMarker extends ConsensusMarker
{
    @Override
    public ConsensusType consensusType(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record)
    {
        int minQual = minQualInMicrosatellite(refGenomeMicrosatellite, record);
        if(minQual > LOW_QUAL_CUTOFF)
            return HIGH_QUAL;

        return NONE;
    }
}
