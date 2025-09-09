package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.minQual;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_SIMPLEX_QUAL;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public final class SbxUtils
{
    public static ConsensusType determineConsensusType(final VariantReadContext readContext, int readVarIndex, final SAMRecord record)
    {
        ConsensusType consensusType = extractConsensusType(record);

        if(consensusType != DUAL)
            return consensusType;

        int duplexBaseIndex = SbxBamUtils.extractDuplexBaseIndex(record);

        int readIndexStart = readVarIndex - readContext.leftCoreLength();
        int readIndexEnd = readVarIndex + readContext.rightCoreLength();

        for(int i = readIndexStart; i <= readIndexEnd; ++i)
        {
            if(!SbxBamUtils.inDuplexRegion(!record.getReadNegativeStrandFlag(), duplexBaseIndex, i))
                return SINGLE;
        }

        return consensusType;
    }

    public static double indelCoreQuality(final VariantReadContext readContext, final SAMRecord record, final int readVarIndex)
    {
        int readIndexStart = max(readVarIndex - readContext.leftCoreLength(), 0);
        int readIndexEnd = min(readVarIndex + readContext.rightCoreLength(), record.getReadBases().length - 1);

        int baseLength = readIndexEnd - readIndexStart + 1;

        if(baseLength <= 0)
            return 0;

        double quality = 0;
        boolean hasMedium = false;

        for(int i = readIndexStart; i <= readIndexEnd; i++)
        {
            quality += record.getBaseQualities()[i];
            hasMedium |= record.getBaseQualities()[i] == SBX_SIMPLEX_QUAL;
        }

        double avgQual = (int)round(quality / baseLength);

        return hasMedium && avgQual > SBX_SIMPLEX_QUAL ? SBX_SIMPLEX_QUAL : avgQual;
    }
}
