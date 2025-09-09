package com.hartwig.hmftools.sage.seqtech;

import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;

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
}
