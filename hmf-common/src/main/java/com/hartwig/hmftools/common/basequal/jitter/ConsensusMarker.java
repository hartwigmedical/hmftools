package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;

import com.hartwig.hmftools.common.sequencing.SequencingType;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public abstract class ConsensusMarker
{
    public abstract ConsensusType consensusType(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record);

    @Nullable
    public static ConsensusMarker fromSequencingType(final SequencingType sequencingType)
    {
        if(sequencingType == SBX)
        {
            return new SBXConsensusMarker();
        }

        return null;
    }
}
