package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.samtools.CigarUtils.cigarBaseLength;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_CONSENSUS_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.DuplicateGroups.calcBaseQualAverage;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_FAIL;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_MISMATCH;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.SUPPLEMENTARY;
import static com.hartwig.hmftools.markdups.umi.UmiConfig.READ_ID_DELIM;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.samtools.CigarUtils;

import htsjdk.samtools.Cigar;
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
        {
            SAMRecord consensusRead = copyPrimaryRead(reads, groupIdentifier);
            return new ConsensusReadInfo(consensusRead, SUPPLEMENTARY);
        }

        boolean isForward = !reads.get(0).getReadNegativeStrandFlag();
        int maxBaseLength = 0;
        boolean hasIndels = false;
        boolean unmapped = false;

        // work out the outermost boundaries - soft-clipped and aligned - from amongst all reads
        ConsensusState consensusState = new ConsensusState(isForward, reads.get(0).getContig());

        for(SAMRecord read : reads)
        {
            maxBaseLength = max(maxBaseLength, read.getReadBases().length);

            if(!read.getReadUnmappedFlag())
                hasIndels |= read.getCigar().getCigarElements().stream().anyMatch(x -> x.getOperator() == I || x.getOperator() == D);
            else
                unmapped = true;

            consensusState.setBoundaries(read);
            consensusState.MapQuality = max(consensusState.MapQuality, read.getMappingQuality());
        }

        if(hasIndels)
        {
            mIndelConsensusReads.buildIndelComponents(reads,  consensusState);

            if(consensusState.outcome() == INDEL_FAIL)
            {
                SAMRecord consensusRead = copyPrimaryRead(reads, groupIdentifier);

                StringJoiner sj = new StringJoiner(", ");
                reads.forEach(x -> sj.add(x.getCigarString()));
                MD_LOGGER.debug("consensus indel mismatch: reads({}) cigars({}) details({})",
                        reads.size(), sj.toString(), readToString(consensusRead));

                return new ConsensusReadInfo(consensusRead, consensusState.outcome());
            }
        }
        else
        {
            consensusState.setBaseLength(maxBaseLength);
            mBaseBuilder.buildReadBases(reads, consensusState);
            consensusState.setOutcome(ALIGNMENT_ONLY);

            if(!unmapped)
                buildCigar(consensusState);
        }

        SAMRecord consensusRead = createConsensusRead(consensusState, reads, groupIdentifier);

        if(mBaseBuilder.config().Debug)
            consensusRead.setMappingQuality(0);

        // sanity check the read
        if(!isValidConsensusRead(consensusRead))
        {
            MD_LOGGER.error("invalid consensus read: umi({}) read: {}", groupIdentifier, readToString(consensusRead));
        }

        return new ConsensusReadInfo(consensusRead, consensusState.outcome());
    }

    private boolean isValidConsensusRead(final SAMRecord consensusRead)
    {
        int baseLength = consensusRead.getReadBases().length;

        if(consensusRead.getBaseQualities().length != baseLength)
            return false;

        if(cigarBaseLength(consensusRead.getCigar()) != baseLength)
            return false;

        return true;
    }

    protected static String formReadId(final String templateReadId, final String groupIdentifier)
    {
        int lastDelim = templateReadId.lastIndexOf(READ_ID_DELIM);
        return lastDelim > 0 ? templateReadId.substring(0, lastDelim) + READ_ID_DELIM + "CNS_" + groupIdentifier
                : templateReadId + READ_ID_DELIM + "CNS_" + groupIdentifier;
    }

    public static SAMRecord createConsensusRead(final ConsensusState state, final List<SAMRecord> reads, final String groupIdentifier)
    {
        SAMRecord initialRead = reads.get(0);
        SAMRecord record = new SAMRecord(initialRead.getHeader());

        record.setReadName(formReadId(initialRead.getReadName(), groupIdentifier));
        record.setReadBases(state.Bases);
        record.setBaseQualities(state.BaseQualities);
        record.setMappingQuality(state.MapQuality);
        record.setReferenceName(initialRead.getReferenceName());

        record.setAlignmentStart(state.MinAlignedPosStart);

        if(!initialRead.getReadUnmappedFlag())
            record.setCigar(new Cigar(state.CigarElements));
        else
            record.setCigar(initialRead.getCigar());

        if(initialRead.getMateReferenceIndex() >= 0)
        {
            record.setMateReferenceName(initialRead.getMateReferenceName());
            record.setMateAlignmentStart(initialRead.getMateAlignmentStart());
            record.setMateReferenceIndex(initialRead.getMateReferenceIndex());
            record.setReadPairedFlag(true);
            record.setProperPairFlag(true);
        }
        else
        {
            record.setReadPairedFlag(false);
            record.setProperPairFlag(false);
        }

        record.setFlags(initialRead.getFlags());
        record.setDuplicateReadFlag(false); // being the new primary
        record.setAttribute(UMI_CONSENSUS_ATTRIBUTE, reads.size());

        if(initialRead.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, initialRead.getAttribute(SUPPLEMENTARY_ATTRIBUTE));
        }

        record.setInferredInsertSize(initialRead.getInferredInsertSize());

        return record;
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
        record.setMappingQuality(primaryRead.getMappingQuality());

        record.setAlignmentStart(primaryRead.getAlignmentStart());
        record.setCigar(primaryRead.getCigar());
        record.setMateReferenceName(primaryRead.getMateReferenceName());
        record.setMateAlignmentStart(primaryRead.getMateAlignmentStart());
        record.setMateReferenceIndex(primaryRead.getMateReferenceIndex());
        record.setFlags(primaryRead.getFlags());
        record.setDuplicateReadFlag(false);

        primaryRead.getAttributes().forEach(x -> record.setAttribute(x.tag, x.value));
        record.setAttribute(UMI_CONSENSUS_ATTRIBUTE, reads.size());

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
