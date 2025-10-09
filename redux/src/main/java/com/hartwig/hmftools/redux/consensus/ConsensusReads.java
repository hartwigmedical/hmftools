package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_READ_INDEX_TAG;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.isSbx;
import static com.hartwig.hmftools.redux.ReduxConstants.CONSENSUS_MAX_DEPTH;
import static com.hartwig.hmftools.redux.ReduxConstants.CONSENSUS_PREFIX;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_FAIL;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.SUPPLEMENTARY;
import static com.hartwig.hmftools.redux.consensus.IndelConsensusReads.selectPrimaryRead;
import static com.hartwig.hmftools.redux.consensus.ReadValidReason.INVALID_BASE;
import static com.hartwig.hmftools.redux.consensus.ReadValidReason.getInvalidBases;
import static com.hartwig.hmftools.redux.consensus.ReadValidReason.isValidRead;
import static com.hartwig.hmftools.redux.duplicate.UmiConfig.READ_ID_DELIM;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.List;
import java.util.Map;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.duplicate.FragmentCoords;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

public class ConsensusReads
{
    private final RefGenome mRefGenome;
    private final BaseBuilder mBaseBuilder;
    private final IndelConsensusReads mIndelConsensusReads;

    private final ConsensusStatistics mConsensusStats;

    public ConsensusReads(final RefGenomeInterface refGenome, final SequencingType sequencingType, final ConsensusStatistics consensusStats)
    {
        mRefGenome = new RefGenome(refGenome);

        mBaseBuilder = new BaseBuilder(mRefGenome, consensusStats, sequencingType);
        mIndelConsensusReads = new IndelConsensusReads(mBaseBuilder);

        mConsensusStats = consensusStats;
    }

    @VisibleForTesting
    public ConsensusReads(final RefGenomeInterface refGenome, final SequencingType sequencingType)
    {
        this(refGenome, sequencingType, new ConsensusStatistics());
    }

    @VisibleForTesting
    public ConsensusReads(final RefGenomeInterface refGenome)
    {
        this(refGenome, ILLUMINA);
    }

    public ConsensusReadInfo createConsensusRead(
            final List<SAMRecord> reads, final FragmentCoords fragmentCoords, @Nullable final String umiId)
    {
        String consensusReadId  = "";

        SAMRecord templateRead = TemplateReads.selectTemplateRead(reads, fragmentCoords);
        consensusReadId = formConsensusReadId(templateRead, umiId);

        if(reads.size() <= 1 || reads.get(0).getReadUnmappedFlag())
        {
            SAMRecord consensusRead = buildFromRead(templateRead, consensusReadId, firstInPair(templateRead));

            return new ConsensusReadInfo(consensusRead, templateRead, SUPPLEMENTARY);
        }

        List<SAMRecord> readsView;

        if(reads.size() < CONSENSUS_MAX_DEPTH)
        {
            readsView = reads;
        }
        else
        {
            readsView = Lists.newArrayList(reads.subList(0, CONSENSUS_MAX_DEPTH));

            if(readsView.stream().noneMatch(x -> x == templateRead)) // ensure it is included since it drives cigar selection
                readsView.add(templateRead);
        }

        boolean isForward = !templateRead.getReadNegativeStrandFlag();
        boolean hasIndels = false;

        // work out the outermost boundaries - clipped and aligned - from amongst all reads
        ConsensusState consensusState = new ConsensusState(consensusReadId, isForward, templateRead.getContig(), mRefGenome);

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

                logInvalidConsensusRead(readsView, null, consensusState, INDEL_FAIL.toString());

                // fall-back to selecting the read with the longest aligned bases, highest average qual
                SAMRecord primaryRead = selectPrimaryRead(readsView);
                SAMRecord consensusRead = buildFromRead(primaryRead, consensusReadId, SamRecordUtils.firstInPair(templateRead));

                return new ConsensusReadInfo(consensusRead, templateRead, consensusState.outcome());
            }
        }
        else
        {
            consensusState.setBaseLength(templateRead.getBaseQualities().length);
            consensusState.setBoundaries(templateRead);
            mBaseBuilder.buildReadBases(readsView, consensusState);
            consensusState.setOutcome(ALIGNMENT_ONLY);

            consensusState.CigarElements.addAll(templateRead.getCigar().getCigarElements());
        }

        mConsensusStats.registerOutcome(consensusState.outcome());

        consensusState.setNumMutations();
        SAMRecord consensusRead = createConsensusRead(consensusState, templateRead);

        if(isSbx())
        {
            int firstDuplexBaseIndex = SbxRoutines.findMaxDuplexBaseIndex(readsView);

            if(firstDuplexBaseIndex >= 0)
                consensusRead.setAttribute(SBX_DUPLEX_READ_INDEX_TAG, firstDuplexBaseIndex);
        }

        if(ReduxConfig.RunChecks)
        {
            ReadValidReason validReason = isValidRead(consensusRead);
            if(validReason != ReadValidReason.OK)
            {
                logInvalidConsensusRead(readsView, consensusRead, consensusState, validReason.toString());
            }
        }

        return new ConsensusReadInfo(consensusRead, templateRead, consensusState.outcome());
    }

    public void setChromosomeLength(int chromosomeLength)
    {
        mBaseBuilder.setChromosomLength(chromosomeLength);
    }

    private static SAMRecord createConsensusRead(final ConsensusState state, final SAMRecord templateRead)
    {
        SAMRecord record = new SAMRecord(templateRead.getHeader());

        record.setReadName(state.ReadId);
        record.setReadBases(state.Bases);
        record.setBaseQualities(state.BaseQualities);
        record.setMappingQuality(state.MapQuality);
        record.setReferenceName(templateRead.getReferenceName());

        record.setAlignmentStart(state.AlignmentStart);

        if(!templateRead.getReadUnmappedFlag())
            record.setCigar(new Cigar(state.CigarElements));
        else
            record.setCigar(templateRead.getCigar());

        templateRead.getAttributes().forEach(x -> record.setAttribute(x.tag, x.value));
        record.setFlags(templateRead.getFlags());
        record.setDuplicateReadFlag(false); // being the new primary
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, state.NumMutations);

        for(Map.Entry<String, Object> tagAndValue : state.Attributes.entrySet())
        {
            String tag = tagAndValue.getKey();
            Object value = tagAndValue.getValue();
            record.setAttribute(tag, value);
        }

        if(!record.getReadPairedFlag())
            return record;

        if(templateRead.getMateReferenceIndex() >= 0)
        {
            record.setMateReferenceName(templateRead.getMateReferenceName());
            record.setMateAlignmentStart(templateRead.getMateAlignmentStart());
            record.setMateReferenceIndex(templateRead.getMateReferenceIndex());
            // mate cigar set with the attributes
        }
        else
        {
            record.setMateAlignmentStart(templateRead.getMateAlignmentStart());
        }

        record.setInferredInsertSize(templateRead.getInferredInsertSize());
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

    public SAMRecord buildFromRead(final SAMRecord read, final String groupReadId, boolean isFirstOfPair)
    {
        SAMRecord record = new SAMRecord(read.getHeader());

        record.setReadName(groupReadId);
        record.setReadBases(read.getReadBases());
        record.setBaseQualities(read.getBaseQualities());
        record.setReferenceName(read.getReferenceName());
        record.setMappingQuality(read.getMappingQuality());

        read.getAttributes().forEach(x -> record.setAttribute(x.tag, x.value));

        record.setAlignmentStart(read.getAlignmentStart());
        record.setCigar(read.getCigar());
        record.setFlags(read.getFlags());
        record.setFirstOfPairFlag(isFirstOfPair);
        record.setSecondOfPairFlag(!isFirstOfPair);
        record.setDuplicateReadFlag(false);
        if(!record.getReadPairedFlag())
            return record;

        record.setMateReferenceName(read.getMateReferenceName());
        record.setMateAlignmentStart(read.getMateAlignmentStart());
        record.setMateReferenceIndex(read.getMateReferenceIndex());

        record.setInferredInsertSize(read.getInferredInsertSize());

        return record;
    }

    private void logInvalidConsensusRead(
            final List<SAMRecord> reads, final SAMRecord consensusRead, final ConsensusState consensusState, final String reason)
    {
        if(!ReduxConfig.RunChecks)
            return;

        String logReason = reason;

        if(reason.equals(INVALID_BASE.toString()))
        {
            List<Integer> invalidBases = getInvalidBases(consensusRead);
            logReason += format(" zero-bases%s", invalidBases);
        }

        RD_LOGGER.error("invalid consensus read({}): readCount({}) {} read: {}",
                logReason, reads.size(), consensusState.IsForward ? "forward" : "reverse",
                consensusRead != null ? readToString(consensusRead) : format("%s:none", consensusState.ReadId));

        // reads
        for(int i = 0; i < reads.size(); ++i)
        {
            RD_LOGGER.debug("read {}: {}", i, readToString(reads.get(i)));
        }
    }
}
