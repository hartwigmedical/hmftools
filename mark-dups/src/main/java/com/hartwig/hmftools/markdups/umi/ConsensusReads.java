package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.markdups.common.DuplicateGroups.calcBaseQualAverage;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_FAIL;
import static com.hartwig.hmftools.markdups.umi.ConsensusState.formReadId;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ConsensusReads
{
    private final BaseBuilder mBaseBuilder;
    private final IndelConsensusReads mIndelConsensusReads;

    public ConsensusReads(final UmiConfig config, final RefGenomeInterface refGenome)
    {
        mBaseBuilder = new BaseBuilder(config, refGenome);
        mIndelConsensusReads = new IndelConsensusReads(mBaseBuilder);
    }

    public ConsensusReadInfo createConsensusRead(final List<SAMRecord> reads, final String groupIdentifier)
    {
        if(reads.size() <= 1)
            return null;

        boolean isForward = !reads.get(0).getReadNegativeStrandFlag();
        int maxBaseLength = 0;
        boolean hasIndels = false;

        // work out the outermost boundaries - soft-clipped and aligned - from amongst all reads
        ConsensusState consensusState = new ConsensusState(isForward);

        for(SAMRecord read : reads)
        {
            maxBaseLength = max(maxBaseLength, read.getReadBases().length);

            hasIndels |= read.getCigar().getCigarElements().stream().anyMatch(x -> x.getOperator() == I || x.getOperator() == D);

            consensusState.setBoundaries(read);
        }

        consensusState.setBaseLength(maxBaseLength);

        if(hasIndels)
        {
            mIndelConsensusReads.buildIndelComponents(reads,  consensusState);

            if(consensusState.outcome() == INDEL_FAIL)
            {
                SAMRecord consensusRead = copyPrimaryRead(reads, groupIdentifier);
                return new ConsensusReadInfo(consensusRead, consensusState.outcome());
            }
        }
        else
        {
            mBaseBuilder.buildReadBases(reads, consensusState);
            consensusState.setOutcome(ALIGNMENT_ONLY);

            buildCigar(consensusState);
        }

        SAMRecord consensusRead = consensusState.createConsensusRead(reads.get(0), groupIdentifier);

        return new ConsensusReadInfo(consensusRead, consensusState.outcome());
    }

    public SAMRecord copyPrimaryRead(final List<SAMRecord> reads, final String groupIdentifier)
    {
        SAMRecord primaryRead = null;

        double maxBaseQual = 0;

        for(SAMRecord read : reads)
        {
            double avgBaseQual = calcBaseQualAverage(read);

            if(avgBaseQual > maxBaseQual)
            {
                maxBaseQual = avgBaseQual;
                primaryRead = read;
            }
        }

        SAMRecord record = new SAMRecord(primaryRead.getHeader());

        record.setReadName(formReadId(primaryRead.getReadName(), groupIdentifier));
        record.setReadBases(primaryRead.getReadBases());
        record.setBaseQualities(primaryRead.getBaseQualities());
        record.setReferenceName(primaryRead.getReferenceName());

        record.setAlignmentStart(primaryRead.getAlignmentStart());
        record.setCigar(primaryRead.getCigar());
        record.setMateReferenceName(primaryRead.getMateReferenceName());
        record.setMateAlignmentStart(primaryRead.getMateAlignmentStart());
        record.setMateReferenceIndex(primaryRead.getMateReferenceIndex());
        record.setFlags(primaryRead.getFlags());
        record.setDuplicateReadFlag(false);

        primaryRead.getAttributes().forEach(x -> record.setAttribute(x.tag, x.value));
        record.setInferredInsertSize(primaryRead.getInferredInsertSize());

        return record;
    }

    private static void buildCigar(final ConsensusState consensusState)
    {
        // build CIGAR from matched and any soft-clipped elements
        int leftSoftClipBases = consensusState.MinAlignedPosStart - consensusState.MinUnclippedPosStart;
        int rightSoftClipBases = consensusState.MaxUnclippedPosEnd - consensusState.MaxAlignedPosEnd;
        int alignedBases = consensusState.MaxAlignedPosEnd - consensusState.MinAlignedPosStart + 1;

        if(leftSoftClipBases > 0)
            consensusState.CigarElements.add(new CigarElement(leftSoftClipBases, CigarOperator.S));

        consensusState.CigarElements.add(new CigarElement(alignedBases, CigarOperator.M));

        if(rightSoftClipBases > 0)
            consensusState.CigarElements.add(new CigarElement(rightSoftClipBases, CigarOperator.S));
    }

}
