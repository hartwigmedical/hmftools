package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.minBaseQualAcrossRange;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public final class SbxUtils
{
    public static final byte MATCHING_BASE_QUALITY_SBX = 10;
    public static final int MQF_NM_1_THRESHOLD_DEDUCTION = 10;

    public static final double GERMLINE_VAF_INDEL_REPEAT_THRESHOLD = 3.5;
    public static final double GERMLINE_VAF_INDEL_REPEAT_MIN_THRESHOLD_FACTOR = 2;
    public static final double GERMLINE_VAF_REL_QUAL_RATIO_THRESHOLD = 0.2;

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

        return minBaseQualAcrossRange(readIndexStart, readIndexEnd, record);
    }
}
