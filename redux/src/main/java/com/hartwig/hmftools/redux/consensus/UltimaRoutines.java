package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;

import htsjdk.samtools.SAMRecord;

public final class UltimaRoutines
{
    public static void finaliseRead(final RefGenomeInterface refGenome, final SAMRecord record)
    {
        ConsensusType consensusType = UltimaBamUtils.deriveConsensusType(record);
        record.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, consensusType.toString());
    }
}
