package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getRequiredSamTag;
import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.emit.BamRecordEmitter;
import com.hartwig.hmftools.tars.liftback.emit.LiftedMatePair;
import com.hartwig.hmftools.tars.liftback.merge.SupplementaryAlignmentResolver;
import com.hartwig.hmftools.tars.liftback.merge.SupplementaryMerger;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGate;
import com.hartwig.hmftools.tars.liftback.select.AlignmentSelector;

import htsjdk.samtools.SAMRecord;

public class LiftBackGroupProcessor
{
    private final AlignmentSelector mAlignmentSelector;
    private final SupplementaryAlignmentResolver mSupplementaryResolver;
    private final OverhangGate mOverhangGate;
    private final GenomicAlignmentScorer mAlignmentScorer;
    private final BamRecordEmitter mEmitter;
    private final LiftBackStats mStats;

    public LiftBackGroupProcessor(
            final AlignmentSelector alignmentSelector, final SupplementaryMerger supplementaryMerger,
            final OverhangGate overhangGate, final GenomicAlignmentScorer alignmentScorer,
            final RefGenomeInterface refGenome, final ExcludedRegions excludedRegions)
    {
        mAlignmentSelector = alignmentSelector;
        mOverhangGate = overhangGate;
        mAlignmentScorer = alignmentScorer;
        mStats = new LiftBackStats();
        mSupplementaryResolver = new SupplementaryAlignmentResolver(
                supplementaryMerger, alignmentScorer, mStats);
        mEmitter = new BamRecordEmitter(refGenome, excludedRegions, mStats);
    }

    public LiftBackStats stats() { return mStats; }

    public void processNameGroup(final List<SAMRecord> group, final Consumer<SAMRecord> consumer)
    {
        mStats.RecordsSeen += group.size();

        List<SAMRecord> firstOfPair = new ArrayList<>();
        List<SAMRecord> secondOfPair = new ArrayList<>();
        for(SAMRecord record : group)
        {
            (firstInPair(record) ? firstOfPair : secondOfPair).add(record);
        }

        LiftedRecord firstAlignments = firstOfPair.isEmpty()
                ? null : mAlignmentSelector.liftPrimaryAlignments(firstOfPair.get(0), mOverhangGate);
        LiftedRecord secondAlignments = secondOfPair.isEmpty()
                ? null : mAlignmentSelector.liftPrimaryAlignments(secondOfPair.get(0), mOverhangGate);

        PreparedRead firstPrepared = prepareRead(firstOfPair, List.of(), firstAlignments);
        List<ChrBaseRegion> firstIntrons = provisionalIntrons(firstPrepared, secondAlignments);
        PreparedRead secondPrepared = prepareRead(secondOfPair, firstIntrons, secondAlignments);

        AlignmentSelector.PairSelection pairSelection = chooseMatePair(firstPrepared, secondPrepared);
        ReadDecision firstDecision = finishRead(
                firstPrepared, pairSelection != null ? pairSelection.first() : null, secondPrepared.primaryAlignments());

        ReadDecision secondDecision = finishRead(
                secondPrepared, pairSelection != null ? pairSelection.second() : null,
                chosenPlacement(firstDecision.primaryResult(), firstAlignments));

        LiftedMatePair matePair = new LiftedMatePair(firstDecision.primaryResult(), secondDecision.primaryResult());
        mEmitter.emit(
                firstOfPair, firstDecision.liftedRecords(), firstDecision.absorbedSupplementaries(),
                firstDecision.primaryUnmapped(), matePair, consumer);
        mEmitter.emit(
                secondOfPair, secondDecision.liftedRecords(), secondDecision.absorbedSupplementaries(),
                secondDecision.primaryUnmapped(), matePair, consumer);
    }

    private PreparedRead prepareRead(
            final List<SAMRecord> records, final List<ChrBaseRegion> mateHintIntrons, final LiftedRecord liftedPrimary)
    {
        if(records.isEmpty())
        {
            return PreparedRead.EMPTY;
        }

        SAMRecord primary = records.get(0);
        if(primary.getSupplementaryAlignmentFlag())
        {
            throw new IllegalStateException(String.format(
                    "read %s: first record in group is a supplementary, expected the primary first",
                    primary.getReadName()));
        }

        List<LiftedRecord> resolved = new ArrayList<>(Collections.nCopies(records.size(), null));
        LiftedRecord primaryAlignments = liftedPrimary;
        boolean hasSupplementaries = records.size() > 1;

        if(hasSupplementaries)
        {
            for(int i = 1; i < records.size(); ++i)
            {
                resolved.set(i, mAlignmentSelector.liftSupplementaryAlignment(records.get(i), mOverhangGate));
            }
            primaryAlignments = mSupplementaryResolver.addSupportedAlignments(
                    records, primaryAlignments, resolved, mateHintIntrons);
        }

        // Pair selection only considers placements that pass the AS floor.
        resolved.set(0, scorePrimaryAlignments(primary, primaryAlignments));
        return new PreparedRead(records, List.copyOf(resolved));
    }

    private ReadDecision finishRead(
            final PreparedRead prepared, final AlignmentSelector.Selection pairSelection, final LiftedRecord mate)
    {
        if(prepared.records().isEmpty())
        {
            return ReadDecision.EMPTY;
        }

        List<SAMRecord> records = prepared.records();
        SAMRecord primary = records.get(0);
        List<LiftedRecord> resolved = new ArrayList<>(prepared.liftedRecords());
        LiftedRecord primaryAlignments = prepared.primaryAlignments();
        if(primaryAlignments.hasPlacement())
        {
            LiftedRecord selected = pairSelection != null
                    ? mAlignmentSelector.selectPrimaryAlignment(primary, primaryAlignments.liftedAlignments(), pairSelection)
                    : mAlignmentSelector.selectPrimaryAlignment(primary, primaryAlignments.liftedAlignments(), mate);
            resolved.set(0, selected);
        }

        SupplementaryAlignmentResolver.Resolution supplementary = mSupplementaryResolver.finish(
                records, List.copyOf(resolved), resolved.get(0));
        List<LiftedRecord> finalRecords = supplementary.liftedRecords();

        UnmapDecision unmapDecision = applyUnmapPolicy(primary, finalRecords.get(0));
        if(unmapDecision.result() != finalRecords.get(0))
        {
            finalRecords = replacePrimary(finalRecords, unmapDecision.result());
        }
        return new ReadDecision(
                finalRecords, supplementary.absorbedSupplementaries(),
                unmapDecision.unmapped());
    }

    private AlignmentSelector.PairSelection chooseMatePair(
            final PreparedRead first, final PreparedRead second)
    {
        if(!first.hasPlacement() || !second.hasPlacement())
        {
            return null;
        }
        return mAlignmentSelector.chooseMatePair(
                first.primary(), first.primaryAlignments().liftedAlignments(),
                second.primary(), second.primaryAlignments().liftedAlignments());
    }

    private List<ChrBaseRegion> provisionalIntrons(final PreparedRead read, final LiftedRecord mate)
    {
        if(!read.hasPlacement())
        {
            return List.of();
        }
        AlignmentSelector.Selection provisional = mAlignmentSelector.choosePrimaryAlignment(
                read.primary(), read.primaryAlignments().liftedAlignments(), mate);
        return provisional.alignment().mergedSupplementaryIntrons();
    }

    private LiftedRecord scorePrimaryAlignments(final SAMRecord primary, final LiftedRecord alignments)
    {
        if(!alignments.hasPlacement()
                || !AlignmentSelector.usesGenomicScore(alignments.liftedAlignments(), primary.getMappingQuality()))
        {
            return alignments;
        }
        return alignments.withLiftedAlignments(mAlignmentScorer.scorePlacements(alignments.liftedAlignments(), primary));
    }

    private LiftedRecord scoreFinalAlignment(final SAMRecord primary, final LiftedRecord lifted)
    {
        if(!lifted.hasPlacement())
        {
            return lifted;
        }

        LiftedAlignment finalAlignment = lifted.primaryAlignment();
        if(finalAlignment.genomicScore() != null || !GenomicAlignmentScorer.requiresRescore(primary, lifted))
        {
            return lifted;
        }

        int score = mAlignmentScorer.scoreRecord(finalAlignment, primary);
        return lifted.withPrimaryAlignment(finalAlignment.withGenomicScore(score));
    }

    private static int finalAlignmentScore(final SAMRecord primary, final LiftedRecord lifted)
    {
        if(lifted.hasPlacement() && lifted.primaryAlignment().genomicScore() != null)
        {
            return lifted.primaryAlignment().genomicScore();
        }

        return getRequiredSamTag(primary, ALIGNMENT_SCORE_ATTRIBUTE, SAMRecord::getIntegerAttribute);
    }

    private UnmapDecision applyUnmapPolicy(final SAMRecord primary, final LiftedRecord result)
    {
        if(mEmitter.excludes(result))
        {
            ++mStats.UnmappedExcludedRegion;
            return new UnmapDecision(LiftedRecord.unmapped(), true);
        }
        if(primary.getReadUnmappedFlag())
        {
            return new UnmapDecision(result, false);
        }
        if(BamRecordEmitter.exceedsMappingCap(primary, result))
        {
            ++mStats.UnmappedOverCap;
            TARS_LOGGER.trace(
                    "over-cap unmap {}:{}: inputMapQuality=0, no XA",
                    primary.getReferenceName(), primary.getAlignmentStart());
            return new UnmapDecision(LiftedRecord.unmapped(), true);
        }

        LiftedRecord scored = scoreFinalAlignment(primary, result);
        int alignmentScore = finalAlignmentScore(primary, scored);
        if(alignmentScore < PRIMARY_AS_UNMAP_THRESHOLD)
        {
            ++mStats.UnmappedLowAlignmentScore;
            TARS_LOGGER.trace(
                    "AS-floor unmap {}: AS={} < {}",
                    primary.getReadName(), alignmentScore, PRIMARY_AS_UNMAP_THRESHOLD);
            return new UnmapDecision(LiftedRecord.unmapped(), true);
        }
        return new UnmapDecision(scored, false);
    }

    private static LiftedRecord chosenPlacement(final LiftedRecord decided, final LiftedRecord fallback)
    {
        if(decided == null || !decided.hasPlacement())
        {
            return fallback;
        }
        return new LiftedRecord(
                decided.updatedMapQuality(), 1, 0, List.of(decided.primaryAlignment()));
    }

    private static List<LiftedRecord> replacePrimary(
            final List<LiftedRecord> records, final LiftedRecord primary)
    {
        List<LiftedRecord> updated = new ArrayList<>(records);
        updated.set(0, primary);
        return updated;
    }

    private record ReadDecision(
            List<LiftedRecord> liftedRecords, Set<Integer> absorbedSupplementaries,
            boolean primaryUnmapped)
    {
        private static final ReadDecision EMPTY = new ReadDecision(List.of(), Set.of(), false);

        private ReadDecision
        {
            liftedRecords = List.copyOf(liftedRecords);
            absorbedSupplementaries = Set.copyOf(absorbedSupplementaries);
        }

        LiftedRecord primaryResult()
        {
            return liftedRecords.isEmpty() ? null : liftedRecords.get(0);
        }
    }

    private record PreparedRead(List<SAMRecord> records, List<LiftedRecord> liftedRecords)
    {
        private static final PreparedRead EMPTY = new PreparedRead(List.of(), List.of());

        private PreparedRead
        {
            records = List.copyOf(records);
            liftedRecords = List.copyOf(liftedRecords);
        }

        SAMRecord primary()
        {
            return records.get(0);
        }

        LiftedRecord primaryAlignments()
        {
            return liftedRecords.isEmpty() ? null : liftedRecords.get(0);
        }

        boolean hasPlacement()
        {
            return primaryAlignments() != null && primaryAlignments().hasPlacement();
        }
    }

    private record UnmapDecision(LiftedRecord result, boolean unmapped)
    {
    }
}
