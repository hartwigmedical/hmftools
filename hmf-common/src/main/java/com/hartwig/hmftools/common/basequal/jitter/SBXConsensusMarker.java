package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.basequal.jitter.ConsensusType.DUPLEX;
import static com.hartwig.hmftools.common.basequal.jitter.ConsensusType.IGNORE;
import static com.hartwig.hmftools.common.basequal.jitter.ConsensusType.NONE;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.SIMPLEX_QUAL;

import htsjdk.samtools.SAMRecord;

public class SBXConsensusMarker extends ConsensusMarker
{
    @Override
    public ConsensusType consensusType(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record)
    {
        int minQual = minQualInMicrosatellite(refGenomeMicrosatellite, record);
        if(minQual == SIMPLEX_QUAL)
            return NONE;

        if(minQual == DUPLEX_QUAL)
            return DUPLEX;

        return IGNORE;
    }
}
