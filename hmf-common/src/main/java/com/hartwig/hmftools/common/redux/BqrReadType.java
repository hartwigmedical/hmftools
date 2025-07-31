package com.hartwig.hmftools.common.redux;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractUmiType;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractConsensusType;

import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.sequencing.UltimaConsensusType;

import htsjdk.samtools.SAMRecord;

public enum BqrReadType
{
    NONE(false),
    SINGLE(false),
    DUAL(true);

    private final boolean mIsHighQuality;

    BqrReadType(boolean isHighQuality) { mIsHighQuality = isHighQuality; }

    public boolean isHighQuality() { return mIsHighQuality; }

    public static BqrReadType extractReadType(final SAMRecord record, final SequencingType sequencingType)
    {
        if(sequencingType == SequencingType.ILLUMINA)
            return BqrReadType.fromUmiType(extractUmiType(record));
        else
            return BqrReadType.fromUltimaType(extractConsensusType(record));
    }

    public static BqrReadType fromUmiType(final UmiReadType umiReadType)
    {
        switch(umiReadType)
        {
            case NONE: return BqrReadType.NONE;
            case SINGLE: return BqrReadType.SINGLE;
            case DUAL: return BqrReadType.DUAL;
            default: return BqrReadType.NONE;
        }
    }

    public static BqrReadType fromUltimaType(final UltimaConsensusType consensusType)
    {
        switch(consensusType)
        {
            case STANDARD: return BqrReadType.NONE;
            case BALANCED: return BqrReadType.DUAL;
            default: return BqrReadType.NONE;
        }
    }
}
