package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;
import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.features.GenomicAlignmentScorer;
import com.hartwig.hmftools.tars.liftback.features.OverhangGate;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryMerger;

import htsjdk.samtools.SAMRecord;

public class LiftBackGroupProcessor
{
    private final PlacementSelector mPlacementSelector;
    private final SupplementaryAlignmentResolver mSupplementaryResolver;
    private final OverhangGate mOverhangGate;
    private final GenomicAlignmentScorer mAlignmentScorer;
    private final BamRecordEmitter mEmitter;
    private final LiftBackStats mStats;

    public LiftBackGroupProcessor(
            final PlacementSelector placementSelector, final SupplementaryMerger supplementaryMerger,
            final OverhangGate overhangGate, final GenomicAlignmentScorer alignmentScorer,
            final RefGenomeInterface refGenome, final ExcludedRegions excludedRegions)
    {
        mPlacementSelector = placementSelector;
        mOverhangGate = overhangGate;
        mAlignmentScorer = alignmentScorer;
        mStats = new LiftBackStats();
        mSupplementaryResolver = new SupplementaryAlignmentResolver(
                supplementaryMerger, alignmentScorer, mStats);
        mEmitter = new BamRecordEmitter(
                placementSelector.contigTranslator(), mSupplementaryResolver.enabled(), refGenome, excludedRegions, mStats);
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

        LiftedRecord firstAlignments = liftPrimaryAlignments(firstOfPair);
        LiftedRecord secondAlignments = liftPrimaryAlignments(secondOfPair);
        LiftedMatePair matePair = new LiftedMatePair();

        PreparedRead firstPrepared = prepareRead(firstOfPair, List.of(), firstAlignments);
        List<ChrBaseRegion> firstIntrons = provisionalIntrons(firstPrepared, secondAlignments);
        PreparedRead secondPrepared = prepareRead(secondOfPair, firstIntrons, secondAlignments);

        PlacementSelector.PairApplyResult pairSelection = chooseMatePair(firstPrepared, secondPrepared);
        ReadDecision firstDecision = finishRead(
                firstPrepared, pairSelection != null ? pairSelection.first() : null, secondPrepared.primaryAlignments());
        recordPrimary(firstOfPair, firstDecision, matePair);

        ReadDecision secondDecision = finishRead(
                secondPrepared, pairSelection != null ? pairSelection.second() : null,
                chosenPlacement(firstDecision.primaryResult(), firstAlignments));
        recordPrimary(secondOfPair, secondDecision, matePair);

        emit(firstOfPair, firstDecision, matePair, consumer);
        emit(secondOfPair, secondDecision, matePair, consumer);
    }

    private LiftedRecord liftPrimaryAlignments(final List<SAMRecord> records)
    {
        if(records.isEmpty() || records.get(0).getReadUnmappedFlag())
        {
            return null;
        }
        return mPlacementSelector.liftPrimaryAlignments(records.get(0), mOverhangGate);
    }

    private PreparedRead prepareRead(
            final List<SAMRecord> records, final List<ChrBaseRegion> mateHintIntrons,
            final LiftedRecord preLiftedPrimaryAlignments)
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
        LiftedRecord primaryAlignments = preLiftedPrimaryAlignments != null
                ? preLiftedPrimaryAlignments : mPlacementSelector.liftPrimaryAlignments(primary, mOverhangGate);
        boolean hasSupplementaries = records.size() > 1;

        if(hasSupplementaries)
        {
            // Step 2: lift the supplementary records and add any spliced primary alignments they support.
            for(int i = 1; i < records.size(); ++i)
            {
                resolved.set(i, mPlacementSelector.liftSupplementaryAlignment(records.get(i), mOverhangGate));
            }
            primaryAlignments = mSupplementaryResolver.addSupportedAlignments(
                    records, primaryAlignments, resolved, mateHintIntrons);
        }

        // Score before pair selection so mate proximity can be applied only to placements that pass the AS floor.
        scorePrimaryAlignments(primary, primaryAlignments);
        resolved.set(0, primaryAlignments);
        return new PreparedRead(records, List.copyOf(resolved));
    }

    private ReadDecision finishRead(
            final PreparedRead prepared, final PlacementSelector.ApplyResult pairSelection, final LiftedRecord mate)
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
                    ? mPlacementSelector.selectPrimaryAlignment(primary, primaryAlignments.liftedAlignments(), pairSelection)
                    : mPlacementSelector.selectPrimaryAlignment(primary, primaryAlignments.liftedAlignments(), mate);
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
                supplementary.introducedIntrons(), unmapDecision.unmapped());
    }

    private PlacementSelector.PairApplyResult chooseMatePair(
            final PreparedRead first, final PreparedRead second)
    {
        if(!first.hasPlacement() || !second.hasPlacement())
        {
            return null;
        }
        return mPlacementSelector.chooseMatePair(
                first.primary(), first.primaryAlignments().liftedAlignments(),
                second.primary(), second.primaryAlignments().liftedAlignments());
    }

    private List<ChrBaseRegion> provisionalIntrons(final PreparedRead read, final LiftedRecord mate)
    {
        if(!read.hasPlacement())
        {
            return List.of();
        }
        PlacementSelector.ApplyResult provisional = mPlacementSelector.choosePrimaryAlignment(
                read.primary(), read.primaryAlignments().liftedAlignments(), mate);
        return provisional.effectivePrimary().MergedSupplementaryIntrons;
    }

    private void scorePrimaryAlignments(final SAMRecord primary, final LiftedRecord alignments)
    {
        if(!alignments.hasPlacement())
        {
            return;
        }
        if(mAlignmentScorer != null
                && PlacementSelector.usesGenomicScore(alignments.liftedAlignments(), primary.getMappingQuality()))
        {
            mAlignmentScorer.scorePlacements(alignments.liftedAlignments(), primary);
        }
    }

    private Integer finalAlignmentScore(final SAMRecord primary, final LiftedRecord lifted)
    {
        // lift-only paths have no supplementary merger and skip the production AS floor
        if(!mSupplementaryResolver.enabled())
        {
            return null;
        }

        if(lifted.hasPlacement())
        {
            LiftedAlignment finalAlignment = lifted.primaryAlignment();

            if(finalAlignment.GenomicScore != Integer.MIN_VALUE)
            {
                return finalAlignment.GenomicScore;
            }

            if(lifted.primaryIndex() != 0 || GenomicAlignmentScorer.requiresRescore(primary, finalAlignment))
            {
                if(mAlignmentScorer != null)
                {
                    int score = mAlignmentScorer.scoreRecord(finalAlignment, primary);
                    if(score != Integer.MIN_VALUE)
                    {
                        finalAlignment.GenomicScore = score;
                        return score;
                    }
                }

                return null;
            }
        }

        return primary.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
    }

    private UnmapDecision applyUnmapPolicy(final SAMRecord primary, final LiftedRecord result)
    {
        if(mEmitter.excludes(result))
        {
            ++mStats.UnmappedExcludedRegion;
            return new UnmapDecision(LiftedRecord.unmapped("excluded_region_unmapped"), true);
        }
        if(primary.getReadUnmappedFlag())
        {
            return new UnmapDecision(result, false);
        }
        if(BamRecordEmitter.exceedsMappingCap(primary, result))
        {
            ++mStats.UnmappedOverCap;
            TARS_LOGGER.trace("over-cap unmap {}:{}: inputMapQuality=0, no XA",
                    primary.getReferenceName(), primary.getAlignmentStart());
            return new UnmapDecision(LiftedRecord.unmapped("over_cap_unmapped"), true);
        }

        Integer alignmentScore = finalAlignmentScore(primary, result);
        if(alignmentScore != null && alignmentScore < PRIMARY_AS_UNMAP_THRESHOLD)
        {
            ++mStats.UnmappedLowAlignmentScore;
            TARS_LOGGER.trace(
                    "AS-floor unmap {}: AS={} < {}",
                    primary.getReadName(), alignmentScore, PRIMARY_AS_UNMAP_THRESHOLD);
            return new UnmapDecision(LiftedRecord.unmapped("low_as_unmapped"), true);
        }
        return new UnmapDecision(result, false);
    }

    private void emit(
            final List<SAMRecord> records, final ReadDecision decision,
            final LiftedMatePair matePair, final Consumer<SAMRecord> consumer)
    {
        mEmitter.emit(
                records, decision.liftedRecords(), decision.absorbedSupplementaries(),
                decision.primaryUnmapped(), matePair, consumer);
    }

    private static void recordPrimary(
            final List<SAMRecord> records, final ReadDecision decision, final LiftedMatePair pair)
    {
        if(!records.isEmpty() && decision.primaryResult() != null)
        {
            pair.recordPrimary(firstInPair(records.get(0)), decision.primaryResult());
        }
    }

    private static LiftedRecord chosenPlacement(final LiftedRecord decided, final LiftedRecord fallback)
    {
        if(decided == null || !decided.hasPlacement())
        {
            return fallback;
        }
        return new LiftedRecord(
                decided.updatedMapQuality(), 1, decided.notes(), 0, List.of(decided.primaryAlignment()));
    }

    private static List<LiftedRecord> replacePrimary(
            final List<LiftedRecord> records, final LiftedRecord primary)
    {
        List<LiftedRecord> updated = new ArrayList<>(records);
        updated.set(0, primary);
        return List.copyOf(updated);
    }

    private record ReadDecision(
            List<LiftedRecord> liftedRecords, Set<Integer> absorbedSupplementaries,
            List<ChrBaseRegion> introducedIntrons, boolean primaryUnmapped)
    {
        private static final ReadDecision EMPTY = new ReadDecision(List.of(), Set.of(), List.of(), false);

        private ReadDecision
        {
            liftedRecords = List.copyOf(liftedRecords);
            absorbedSupplementaries = Set.copyOf(absorbedSupplementaries);
            introducedIntrons = List.copyOf(introducedIntrons);
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
