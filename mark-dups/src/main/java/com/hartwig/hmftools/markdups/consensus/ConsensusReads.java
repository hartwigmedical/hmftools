package com.hartwig.hmftools.markdups.consensus;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.samtools.CigarUtils.cigarBaseLength;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.Constants.CONSENSUS_MAX_DEPTH;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;
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
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ConsensusReads
{
    private final RefGenomeInterface mRefGenome;
    private final BaseBuilder mBaseBuilder;
    private final IndelConsensusReads mIndelConsensusReads;

    private final ConsensusStatistics mConsensusStats;
    private boolean mValidateConsensusReads;

    private static final String CONSENSUS_PREFIX = "CNS_";

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

    public ConsensusReadInfo createConsensusRead(final List<SAMRecord> reads, final String groupIdentifier)
    {
        if(reads.size() <= 1 || reads.get(0).getReadUnmappedFlag())
        {
            SAMRecord consensusRead = copyPrimaryRead(reads.get(0), groupIdentifier);
            return new ConsensusReadInfo(consensusRead, SUPPLEMENTARY);
        }

        List<SAMRecord> readsView = reads.subList(0, min(CONSENSUS_MAX_DEPTH, reads.size()));

        boolean isForward = !readsView.get(0).getReadNegativeStrandFlag();
        boolean hasIndels = false;

        // work out the outermost boundaries - soft-clipped and aligned - from amongst all reads
        ConsensusState consensusState = new ConsensusState(isForward, readsView.get(0).getContig(), mRefGenome);

        for(SAMRecord read : readsView)
        {
            hasIndels |= read.getCigar().getCigarElements().stream().anyMatch(x -> x.getOperator() == I || x.getOperator() == D);
            consensusState.MapQuality = max(consensusState.MapQuality, read.getMappingQuality());
        }

        if(hasIndels)
        {
            mIndelConsensusReads.buildIndelComponents(readsView, consensusState);

            if(consensusState.outcome() == INDEL_FAIL)
            {
                mConsensusStats.registerOutcome(INDEL_FAIL);

                logInvalidConsensusRead(readsView, null, groupIdentifier, consensusState, INDEL_FAIL.toString());

                // fall-back to selecting the read with the longest aligned bases, highest average qual
                SAMRecord primaryRead = selectPrimaryRead(readsView);
                SAMRecord consensusRead = copyPrimaryRead(primaryRead, groupIdentifier);

                return new ConsensusReadInfo(consensusRead, consensusState.outcome());
            }
        }
        else
        {
            Map<String, CigarFrequency> cigarFrequencies = CigarFrequency.buildFrequencies(readsView);
            SAMRecord selectedConsensusRead = cigarFrequencies.size() > 1 ? selectConsensusRead(cigarFrequencies) : readsView.get(0);

            consensusState.setBaseLength(selectedConsensusRead.getBaseQualities().length);
            consensusState.setBoundaries(selectedConsensusRead);
            mBaseBuilder.buildReadBases(readsView, consensusState);
            consensusState.setOutcome(ALIGNMENT_ONLY);

            consensusState.CigarElements.addAll(selectedConsensusRead.getCigar().getCigarElements());
        }

        mConsensusStats.registerOutcome(consensusState.outcome());

        consensusState.setNumMutations();
        SAMRecord consensusRead = createConsensusRead(consensusState, readsView, groupIdentifier);

        if(mValidateConsensusReads)
        {
            ValidationReason validReason = isValidConsensusRead(consensusRead);
            if(validReason != ValidationReason.OK)
            {
                logInvalidConsensusRead(readsView, consensusRead, groupIdentifier, consensusState, validReason.toString());
            }
        }

        return new ConsensusReadInfo(consensusRead, consensusState.outcome());
    }

    protected static SAMRecord selectConsensusRead(final Map<String,CigarFrequency> cigarFrequencies)
    {
        int maxCigarFreq = cigarFrequencies.values().stream().mapToInt(x -> x.Frequency).max().orElse(0);
        List<CigarFrequency> maxCigarFrequencies = cigarFrequencies.values().stream().filter(x -> x.Frequency == maxCigarFreq).collect(Collectors.toList());

        if(maxCigarFrequencies.size() == 1)
            return maxCigarFrequencies.get(0).SampleRead;

        // find the most common read by CIGAR, and where there are equal counts choose the one with the least soft-clips
        SAMRecord selectedRead = null;
        int minScBases = 0;

        for(CigarFrequency cigarFrequency : maxCigarFrequencies)
        {
            int scBases = cigarElementLength(cigarFrequency.SampleRead, S);

            if(selectedRead == null || scBases < minScBases)
            {
                selectedRead = cigarFrequency.SampleRead;
                minScBases = scBases;
            }
        }

        return selectedRead;
    }

    private static int cigarElementLength(final SAMRecord read, final CigarOperator operator)
    {
        return read.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == operator).mapToInt(x -> x.getLength()).sum();
    }

    public void setChromosomeLength(int chromosomeLength)
    {
        mBaseBuilder.setChromosomLength(chromosomeLength);
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

    protected static String formReadId(final String templateReadId, final String groupIdentifier)
    {
        int lastDelim = templateReadId.lastIndexOf(READ_ID_DELIM);
        return lastDelim > 0 ? templateReadId.substring(0, lastDelim) + READ_ID_DELIM + CONSENSUS_PREFIX + groupIdentifier
                : templateReadId + READ_ID_DELIM + CONSENSUS_PREFIX + groupIdentifier;
    }

    private static SAMRecord createConsensusRead(final ConsensusState state, final List<SAMRecord> reads, final String groupIdentifier)
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

        Map<String,CigarFrequency> mateCigarFrequencies = CigarFrequency.buildMateFrequencies(reads);

        SAMRecord selectedMateRead;

        if(mateCigarFrequencies.isEmpty())
            selectedMateRead = initialRead;
        else if(mateCigarFrequencies.size() == 1)
            selectedMateRead = mateCigarFrequencies.values().iterator().next().SampleRead;
        else
            selectedMateRead = selectConsensusRead(mateCigarFrequencies);

        initialRead.getAttributes().forEach(x -> record.setAttribute(x.tag, x.value));

        if(initialRead.getMateReferenceIndex() >= 0)
        {
            record.setMateReferenceName(selectedMateRead.getMateReferenceName());
            record.setMateAlignmentStart(selectedMateRead.getMateAlignmentStart());
            record.setMateReferenceIndex(selectedMateRead.getMateReferenceIndex());
            record.setAttribute(MATE_CIGAR_ATTRIBUTE, selectedMateRead.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
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

        record.setInferredInsertSize(initialRead.getInferredInsertSize());
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, state.NumMutations);
        return record;
    }

    public SAMRecord copyPrimaryRead(final SAMRecord read, final String groupIdentifier)
    {
        SAMRecord record = new SAMRecord(read.getHeader());

        record.setReadName(formReadId(read.getReadName(), groupIdentifier));
        record.setReadBases(read.getReadBases());
        record.setBaseQualities(read.getBaseQualities());
        record.setReferenceName(read.getReferenceName());
        record.setMappingQuality(read.getMappingQuality());

        record.setAlignmentStart(read.getAlignmentStart());
        record.setCigar(read.getCigar());
        record.setMateReferenceName(read.getMateReferenceName());
        record.setMateAlignmentStart(read.getMateAlignmentStart());
        record.setMateReferenceIndex(read.getMateReferenceIndex());
        record.setFlags(read.getFlags());
        record.setDuplicateReadFlag(false);

        read.getAttributes().forEach(x -> record.setAttribute(x.tag, x.value));
        record.setInferredInsertSize(read.getInferredInsertSize());

        return record;
    }
}
