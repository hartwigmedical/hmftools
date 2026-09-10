package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.features.GenomicAlignmentScorer;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryMerger;

import htsjdk.samtools.SAMRecord;

// Owns Step 2 from lifted supplementary placements through merge candidates and final supplementary records.
final class SupplementaryAlignmentResolver
{
    private final SupplementaryMerger mMerger;
    private final GenomicAlignmentScorer mAlignmentScorer;
    private final LiftBackStats mStats;

    SupplementaryAlignmentResolver(
            final SupplementaryMerger merger, final GenomicAlignmentScorer alignmentScorer,
            final LiftBackStats stats)
    {
        mMerger = merger;
        mAlignmentScorer = alignmentScorer;
        mStats = stats;
    }

    boolean enabled()
    {
        return mMerger != null;
    }

    LiftedRecord addSupportedAlignments(
            final List<SAMRecord> records, final LiftedRecord primaryAlignments,
            final List<LiftedRecord> liftedRecords, final List<ChrBaseRegion> mateHintIntrons)
    {
        SAMRecord primary = records.get(0);
        if(!enabled() || primary.getReadUnmappedFlag() || !primaryAlignments.hasPlacement())
        {
            return primaryAlignments;
        }

        List<SupplementaryMerger.Supplementary> supplementaries = placements(records, liftedRecords).alignments();
        if(supplementaries.isEmpty())
        {
            return primaryAlignments;
        }
        ++mStats.MergeableSupplementaries;

        List<LiftedAlignment> expanded = new ArrayList<>(primaryAlignments.liftedAlignments());
        Set<SupportedAlignmentKey> seen = new HashSet<>();
        for(LiftedAlignment alignment : primaryAlignments.liftedAlignments())
        {
            seen.add(new SupportedAlignmentKey(alignment.key(), List.of()));
        }

        for(LiftedAlignment alignment : primaryAlignments.liftedAlignments())
        {
            if(alignment.Dropped || alignment.LiftedCigar == null)
            {
                continue;
            }

            SupplementaryMerger.Result result = mMerger.merge(new SupplementaryMerger.Placement(
                    alignment.LiftedChromosome, alignment.ForwardStrand, primary.getReadLength(),
                    alignment.LiftedPos, alignment.LiftedCigar, supplementaries,
                    primary.getReadBases(), mateHintIntrons));
            ++mStats.SupplementaryMergeAttempts;
            if(!result.merged())
            {
                mStats.recordSupplementaryMergeRejection(result.rejectReason());
                continue;
            }
            ++mStats.SupplementaryMergeSuccessfulCandidates;

            int supplementaryMapQuality = 0;
            for(Integer recordIndex : result.droppedSupplementaryIndices())
            {
                supplementaryMapQuality = Math.max(
                        supplementaryMapQuality, records.get(recordIndex).getMappingQuality());
            }

            LiftedAlignment merged = alignment.withSupplementaryMerge(
                    result.mergedStart(), result.mergedCigar(), result.droppedSupplementaryIndices(),
                    result.introducedIntrons(), supplementaryMapQuality, result.spliceStrand());
            if(seen.add(new SupportedAlignmentKey(merged.key(), result.droppedSupplementaryIndices())))
            {
                expanded.add(merged);
            }
        }

        return primaryAlignments.withLiftedAlignments(List.copyOf(expanded));
    }

    Resolution finish(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords,
            final LiftedRecord primaryResult)
    {
        Set<Integer> absorbed = absorbedSupplementaries(records.size(), primaryResult);
        if(hasSupplementaryMerge(primaryResult))
        {
            ++mStats.SupplementaryMerges;
            mStats.SupplementariesAbsorbed += absorbed.size();
        }

        List<LiftedRecord> selected = selectRecordAlignments(
                records, liftedRecords, primaryResult, absorbed);
        return new Resolution(annotateSpliceStrands(selected), absorbed);
    }

    private List<LiftedRecord> selectRecordAlignments(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords,
            final LiftedRecord primaryResult, final Set<Integer> absorbedSupplementaries)
    {
        if(!enabled() || !primaryResult.hasPlacement())
        {
            return liftedRecords;
        }

        SupplementaryPlacements placements = placements(records, liftedRecords);
        if(placements.alignments().isEmpty())
        {
            return liftedRecords;
        }

        LiftedAlignment primary = primaryResult.primaryAlignment();
        SupplementaryMerger.Placement context = new SupplementaryMerger.Placement(
                primary.LiftedChromosome, primary.ForwardStrand, records.get(0).getReadLength(),
                primary.LiftedPos, primary.LiftedCigar, placements.alignments(), records.get(0).getReadBases(), List.of());
        List<SupplementaryMerger.Supplementary> selected =
                SupplementaryMerger.selectSupplementaryPlacements(context);

        List<LiftedRecord> revised = null;
        for(SupplementaryMerger.Supplementary supplementary : selected)
        {
            int recordIndex = supplementary.index();
            if(absorbedSupplementaries.contains(recordIndex))
            {
                continue;
            }

            Integer alignmentIndex = placements.alignmentIndices().get(supplementary);
            LiftedRecord lifted = liftedRecords.get(recordIndex);
            if(alignmentIndex == null || lifted == null || alignmentIndex == lifted.primaryIndex())
            {
                continue;
            }

            LiftedRecord selectedRecord = lifted.withPrimaryIndex(alignmentIndex);
            int score = mAlignmentScorer != null
                    ? mAlignmentScorer.scoreRecord(selectedRecord.primaryAlignment(), records.get(recordIndex))
                    : Integer.MIN_VALUE;
            if(score != Integer.MIN_VALUE)
            {
                selectedRecord.primaryAlignment().GenomicScore = score;
            }

            if(revised == null)
            {
                revised = new ArrayList<>(liftedRecords);
            }
            revised.set(recordIndex, selectedRecord);
        }

        return revised != null ? List.copyOf(revised) : liftedRecords;
    }

    private static SupplementaryPlacements placements(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords)
    {
        List<SupplementaryMerger.Supplementary> alignments = new ArrayList<>();
        IdentityHashMap<SupplementaryMerger.Supplementary, Integer> alignmentIndices = new IdentityHashMap<>();
        for(int i = 1; i < records.size(); ++i)
        {
            SAMRecord record = records.get(i);
            LiftedRecord lifted = liftedRecords.get(i);
            if(!record.getSupplementaryAlignmentFlag() || record.getReadUnmappedFlag()
                    || lifted == null || !lifted.hasPlacement())
            {
                continue;
            }

            for(int alignmentIndex = 0; alignmentIndex < lifted.liftedAlignments().size(); ++alignmentIndex)
            {
                LiftedAlignment alignment = lifted.liftedAlignments().get(alignmentIndex);
                if(alignment.Dropped || alignment.LiftedCigar == null)
                {
                    continue;
                }

                SupplementaryMerger.Supplementary supplementary = new SupplementaryMerger.Supplementary(
                        i, alignment.LiftedChromosome, alignment.ForwardStrand,
                        alignment.LiftedPos, alignment.LiftedCigar, record.getMappingQuality());
                alignments.add(supplementary);
                alignmentIndices.put(supplementary, alignmentIndex);
            }
        }
        return new SupplementaryPlacements(List.copyOf(alignments), alignmentIndices);
    }

    private List<LiftedRecord> annotateSpliceStrands(final List<LiftedRecord> liftedRecords)
    {
        if(!enabled())
        {
            return liftedRecords;
        }

        List<LiftedRecord> annotated = null;
        for(int i = 0; i < liftedRecords.size(); ++i)
        {
            LiftedRecord lifted = liftedRecords.get(i);
            if(lifted == null || !lifted.hasPlacement() || !lifted.hasNCigar() || lifted.transcriptStrand() != 0)
            {
                continue;
            }

            int strand = mMerger.spliceStrand(
                    lifted.finalChromosome(), lifted.finalPos(), lifted.finalCigar());
            if(strand != 0)
            {
                if(annotated == null)
                {
                    annotated = new ArrayList<>(liftedRecords);
                }
                annotated.set(i, lifted.withPrimaryTranscriptStrand(strand));
            }
        }
        return annotated != null ? List.copyOf(annotated) : liftedRecords;
    }

    private static Set<Integer> absorbedSupplementaries(
            final int recordCount, final LiftedRecord primaryResult)
    {
        if(!hasSupplementaryMerge(primaryResult))
        {
            return Set.of();
        }

        Set<Integer> absorbed = new HashSet<>();
        for(Integer index : primaryResult.primaryAlignment().MergedSupplementaryIndices)
        {
            if(index >= 0 && index < recordCount)
            {
                absorbed.add(index);
            }
        }
        return Set.copyOf(absorbed);
    }

    private static boolean hasSupplementaryMerge(final LiftedRecord primaryResult)
    {
        return primaryResult != null && primaryResult.hasPlacement()
                && primaryResult.primaryAlignment().hasSupplementaryMerge();
    }

    record Resolution(
            List<LiftedRecord> liftedRecords, Set<Integer> absorbedSupplementaries)
    {
        Resolution
        {
            liftedRecords = List.copyOf(liftedRecords);
            absorbedSupplementaries = Set.copyOf(absorbedSupplementaries);
        }
    }

    private record SupportedAlignmentKey(
            AlignmentKey alignment, List<Integer> absorbedSupplementaries)
    {
        private SupportedAlignmentKey
        {
            absorbedSupplementaries = List.copyOf(absorbedSupplementaries);
        }
    }

    private record SupplementaryPlacements(
            List<SupplementaryMerger.Supplementary> alignments,
            IdentityHashMap<SupplementaryMerger.Supplementary, Integer> alignmentIndices)
    {
    }
}
