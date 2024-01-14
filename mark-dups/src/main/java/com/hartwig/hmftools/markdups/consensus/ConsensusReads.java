package com.hartwig.hmftools.markdups.consensus;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.samtools.CigarUtils.cigarBaseLength;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.Constants.CONSENSUS_MAX_DEPTH;
import static com.hartwig.hmftools.markdups.common.Constants.CONSENSUS_PREFIX;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;
import static com.hartwig.hmftools.markdups.consensus.CigarFrequency.selectTemplateRead;
import static com.hartwig.hmftools.markdups.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.markdups.consensus.ConsensusOutcome.INDEL_FAIL;
import static com.hartwig.hmftools.markdups.consensus.ConsensusOutcome.SUPPLEMENTARY;
import static com.hartwig.hmftools.markdups.consensus.IndelConsensusReads.alignedOrSoftClip;
import static com.hartwig.hmftools.markdups.consensus.IndelConsensusReads.selectPrimaryRead;
import static com.hartwig.hmftools.markdups.umi.UmiConfig.READ_ID_DELIM;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ConsensusReads
{
    private final RefGenomeInterface mRefGenome;
    private final BaseBuilder mBaseBuilder;
    private final IndelConsensusReads mIndelConsensusReads;

    private final ConsensusStatistics mConsensusStats;
    private boolean mValidateConsensusReads;

    public ConsensusReads(final RefGenomeInterface refGenome, final ConsensusStatistics consensusStats)
    {
        mRefGenome = refGenome;
        mBaseBuilder = new BaseBuilder(refGenome, consensusStats);
        mConsensusStats = consensusStats;
        mIndelConsensusReads = new IndelConsensusReads(mBaseBuilder);
        mValidateConsensusReads = false;
    }

    @VisibleForTesting
    public ConsensusReads(final RefGenomeInterface refGenome)
    {
        this(refGenome, new ConsensusStatistics());
    }

    public void setDebugOptions(boolean validateConsensusReads)
    {
        mValidateConsensusReads = validateConsensusReads;
    }
    public ConsensusStatistics consensusStats() { return mConsensusStats; }

    public ConsensusReadInfo createConsensusRead(
            final List<SAMRecord> reads, @Nullable final SAMRecord previousTemplateRead,
            @Nullable final String groupReadId, @Nullable final String umiId)
    {
        String consensusReadId;
        SAMRecord templateRead;

        if(previousTemplateRead == null)
        {
            templateRead = selectTemplateRead(reads);
            consensusReadId = formConsensusReadId(templateRead, umiId);
        }
        else
        {
            // match the mate or supplmentary template read to that of the primary
            templateRead = reads.stream().filter(x -> x.getReadName().equals(previousTemplateRead.getReadName())).findFirst().orElse(null);
            consensusReadId = groupReadId;
        }

        if(reads.size() <= 1 || reads.get(0).getReadUnmappedFlag())
        {
            SAMRecord consensusRead = buildFromRead(
                    templateRead != null ? templateRead : reads.get(0),
                    consensusReadId,
                    templateRead == null ? previousTemplateRead : null);

            return new ConsensusReadInfo(consensusRead, templateRead, SUPPLEMENTARY);
        }

        List<SAMRecord> readsView = reads.size() < CONSENSUS_MAX_DEPTH ? reads : reads.subList(0, CONSENSUS_MAX_DEPTH);

        boolean isForward = !templateRead.getReadNegativeStrandFlag();
        boolean hasIndels = false;

        // work out the outermost boundaries - soft-clipped and aligned - from amongst all reads
        ConsensusState consensusState = new ConsensusState(isForward, templateRead.getContig(), mRefGenome);

        for(SAMRecord read : readsView)
        {
            hasIndels |= read.getCigar().getCigarElements().stream().anyMatch(x -> x.getOperator() == I || x.getOperator() == D);
            consensusState.MapQuality = max(consensusState.MapQuality, read.getMappingQuality());
        }

        if(hasIndels)
        {
            mIndelConsensusReads.buildIndelComponents(readsView, consensusState, templateRead);

            if(consensusState.outcome() == INDEL_FAIL)
            {
                mConsensusStats.registerOutcome(INDEL_FAIL);

                logInvalidConsensusRead(readsView, null, consensusReadId, consensusState, INDEL_FAIL.toString());

                // fall-back to selecting the read with the longest aligned bases, highest average qual
                SAMRecord primaryRead = selectPrimaryRead(readsView);
                SAMRecord consensusRead = buildFromRead(primaryRead, consensusReadId, null);

                return new ConsensusReadInfo(consensusRead, templateRead, consensusState.outcome());
            }
        }
        else
        {
            // Map<String, CigarFrequency> cigarFrequencies = CigarFrequency.buildFrequencies(readsView);
            // SAMRecord selectedConsensusRead = cigarFrequencies.size() > 1 ? selectConsensusRead(cigarFrequencies) : readsView.get(0);

            consensusState.setBaseLength(templateRead.getBaseQualities().length);
            consensusState.setBoundaries(templateRead);
            mBaseBuilder.buildReadBases(readsView, consensusState);
            consensusState.setOutcome(ALIGNMENT_ONLY);

            consensusState.CigarElements.addAll(templateRead.getCigar().getCigarElements());
        }

        mConsensusStats.registerOutcome(consensusState.outcome());

        consensusState.setNumMutations();
        SAMRecord consensusRead = createConsensusRead(consensusState, templateRead, consensusReadId);

        if(mValidateConsensusReads)
        {
            ValidationReason validReason = isValidConsensusRead(consensusRead);
            if(validReason != ValidationReason.OK)
            {
                logInvalidConsensusRead(readsView, consensusRead, consensusReadId, consensusState, validReason.toString());
            }
        }

        return new ConsensusReadInfo(consensusRead, templateRead, consensusState.outcome());
    }

    public void setChromosomeLength(int chromosomeLength)
    {
        mBaseBuilder.setChromosomLength(chromosomeLength);
    }

    private static SAMRecord createConsensusRead(final ConsensusState state, final SAMRecord templateRead, final String groupReadId)
    {
        SAMRecord record = new SAMRecord(templateRead.getHeader());

        record.setReadName(groupReadId);
        record.setReadBases(state.Bases);
        record.setBaseQualities(state.BaseQualities);
        record.setMappingQuality(state.MapQuality);
        record.setReferenceName(templateRead.getReferenceName());

        record.setAlignmentStart(state.MinAlignedPosStart);

        if(!templateRead.getReadUnmappedFlag())
            record.setCigar(new Cigar(state.CigarElements));
        else
            record.setCigar(templateRead.getCigar());

        templateRead.getAttributes().forEach(x -> record.setAttribute(x.tag, x.value));

        if(templateRead.getMateReferenceIndex() >= 0)
        {
            record.setMateReferenceName(templateRead.getMateReferenceName());
            record.setMateAlignmentStart(templateRead.getMateAlignmentStart());
            record.setMateReferenceIndex(templateRead.getMateReferenceIndex());
            record.setReadPairedFlag(true);
            record.setProperPairFlag(true);
        }
        else
        {
            record.setReadPairedFlag(false);
            record.setProperPairFlag(false);
        }

        record.setFlags(templateRead.getFlags());
        record.setDuplicateReadFlag(false); // being the new primary

        record.setInferredInsertSize(templateRead.getInferredInsertSize());
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, state.NumMutations);
        return record;
    }

    @VisibleForTesting
    public static String formConsensusReadId(final SAMRecord templateRead, @Nullable final String umiId)
    {
        // take the first read's ID after sorting, include the CNS identifier, and append the UMI if it has one
        String readId = templateRead.getReadName();

        int lastDelim = readId.lastIndexOf(READ_ID_DELIM);

        if(lastDelim <= 0)
        {
            return umiId != null ? readId + READ_ID_DELIM + CONSENSUS_PREFIX + umiId : CONSENSUS_PREFIX + readId;
        }

        String groupId = readId.substring(0, lastDelim) + READ_ID_DELIM + CONSENSUS_PREFIX;

        if(umiId != null)
            return groupId + umiId;
        else
            return groupId + readId.substring(lastDelim + 1);
    }

    public SAMRecord buildFromRead(final SAMRecord read, final String groupReadId, @Nullable final SAMRecord primaryTemplateRead)
    {
        SAMRecord record = new SAMRecord(read.getHeader());

        record.setReadName(groupReadId);
        record.setReadBases(read.getReadBases());
        record.setBaseQualities(read.getBaseQualities());
        record.setReferenceName(read.getReferenceName());
        record.setMappingQuality(read.getMappingQuality());

        read.getAttributes().forEach(x -> record.setAttribute(x.tag, x.value));

        if(read.getReadUnmappedFlag() && primaryTemplateRead != null)
        {
            // rather than use this duplicate's unmapped fields, infer from the mapped primary used for the consensus for consistency
            record.setReadUnmappedFlag(true);
            record.setMateReferenceName(primaryTemplateRead.getReferenceName());
            record.setMateAlignmentStart(primaryTemplateRead.getAlignmentStart());
            record.setMateReferenceIndex(primaryTemplateRead.getReferenceIndex());

            record.setAlignmentStart(primaryTemplateRead.getAlignmentStart());
            record.setCigarString(NO_CIGAR);

            //record.setFlags(read.getFlags());
            record.setAttribute(MATE_CIGAR_ATTRIBUTE, primaryTemplateRead.getCigarString());
            record.setReadPairedFlag(true);

            if(primaryTemplateRead.getFirstOfPairFlag())
                record.setSecondOfPairFlag(true);
            else
                record.setFirstOfPairFlag(true);

            record.setMateNegativeStrandFlag(primaryTemplateRead.getReadNegativeStrandFlag());
            record.setReadNegativeStrandFlag(primaryTemplateRead.getMateNegativeStrandFlag());
        }
        else
        {
            record.setAlignmentStart(read.getAlignmentStart());
            record.setCigar(read.getCigar());
            record.setMateReferenceName(read.getMateReferenceName());
            record.setMateAlignmentStart(read.getMateAlignmentStart());
            record.setMateReferenceIndex(read.getMateReferenceIndex());
            record.setFlags(read.getFlags());
        }

        record.setDuplicateReadFlag(false);

        record.setInferredInsertSize(read.getInferredInsertSize());

        return record;
    }

    private enum ValidationReason
    {
        OK,
        BASE_LENGTH,
        CIGAR_LENGTH,
        CIGAR_ELEMENTS;
    }

    private void logInvalidConsensusRead(
            final List<SAMRecord> reads, final SAMRecord consensusRead, final String groupIdentifier, final ConsensusState consensusState,
            final String reason)
    {
        if(!mValidateConsensusReads)
            return;

        MD_LOGGER.error("invalid consensus read({}): groupId({}) readCount({}) {} read: {}",
                reason, groupIdentifier, reads.size(), consensusState.IsForward ? "forward" : "reverse",
                consensusRead != null ? readToString(consensusRead) : "none");

        // reads
        for(int i = 0; i < reads.size(); ++i)
        {
            MD_LOGGER.debug("read {}: {}", i, readToString(reads.get(i)));
        }
    }

    private ValidationReason isValidConsensusRead(final SAMRecord consensusRead)
    {
        int baseLength = consensusRead.getReadBases().length;

        if(consensusRead.getBaseQualities().length != baseLength)
            return ValidationReason.BASE_LENGTH;

        if(cigarBaseLength(consensusRead.getCigar()) != baseLength)
            return ValidationReason.CIGAR_LENGTH;

        int cigarCount = consensusRead.getCigar().getCigarElements().size();
        if(cigarCount > 1)
        {
            for(int i = 0; i < consensusRead.getCigar().getCigarElements().size(); ++i)
            {
                CigarElement element = consensusRead.getCigar().getCigarElements().get(i);

                if(i == 0 || i == cigarCount - 1)
                {
                    if(!alignedOrSoftClip(element.getOperator()) && element.getOperator() != I)
                        return ValidationReason.CIGAR_ELEMENTS;
                }
                else if(element.getOperator() == S)
                {
                    return ValidationReason.CIGAR_ELEMENTS;
                }
            }
        }
        else if(consensusRead.getCigar().getCigarElements().get(0).getOperator() != M)
        {
            return ValidationReason.CIGAR_ELEMENTS;
        }

        return ValidationReason.OK;
    }
}
