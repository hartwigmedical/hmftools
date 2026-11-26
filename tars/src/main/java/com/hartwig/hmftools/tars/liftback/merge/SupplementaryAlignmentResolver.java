package com.hartwig.hmftools.tars.liftback.merge;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.AlignmentKey;
import com.hartwig.hmftools.tars.liftback.GenomicAlignmentScorer;
import com.hartwig.hmftools.tars.liftback.LiftBackStats;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;
import com.hartwig.hmftools.tars.liftback.LiftedRecord;
import com.hartwig.hmftools.tars.liftback.select.AlignmentSelector;

import htsjdk.samtools.SAMRecord;

// Owns Step 2 from lifted supplementary placements through merge candidates and final supplementary records.
public final class SupplementaryAlignmentResolver
{
    private final SupplementaryMerger mMerger;
    private final GenomicAlignmentScorer mAlignmentScorer;
    private final LiftBackStats mStats;

    public SupplementaryAlignmentResolver(
            final SupplementaryMerger merger, final GenomicAlignmentScorer alignmentScorer,
            final LiftBackStats stats)
    {
        mMerger = merger;
        mAlignmentScorer = alignmentScorer;
        mStats = stats;
    }

    public LiftedRecord addSupportedAlignments(
            final List<SAMRecord> records, final LiftedRecord primaryAlignments,
            final List<LiftedRecord> liftedRecords, final List<ChrBaseRegion> mateHintIntrons)
    {
        SAMRecord primary = records.get(0);
        if(!primaryAlignments.hasPlacement())
        {
            return primaryAlignments;
        }

        List<AlignmentSelector.RecordAlignment> supplementaries = placements(records, liftedRecords);
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
            if(alignment.dropped())
            {
                continue;
            }

            SupplementaryMerger.Result result = mMerger.merge(new SupplementaryMerger.Placement(
                    alignment.liftedPosition().Chromosome, alignment.forwardStrand(), primary.getReadLength(),
                    alignment.liftedPosition().Position, alignment.liftedCigar(), supplementaries,
                    primary.getReadBases(), mateHintIntrons));
            ++mStats.SupplementaryMergeAttempts;
            if(!result.merged())
            {
                mStats.recordSupplementaryMergeRejection(result.rejectReason());
                continue;
            }
            ++mStats.SupplementaryMergeSuccessfulCandidates;

            int supplementaryMapQuality = 0;
            for(Integer recordIndex : result.absorbedSupplementaryIndices())
            {
                supplementaryMapQuality = Math.max(
                        supplementaryMapQuality, records.get(recordIndex).getMappingQuality());
            }

            LiftedAlignment merged = alignment.withSupplementaryMerge(
                    result.mergedStart(), result.mergedCigar(), result.absorbedSupplementaryIndices(),
                    result.introducedIntrons(), supplementaryMapQuality, result.spliceStrand());
            if(seen.add(new SupportedAlignmentKey(merged.key(), result.absorbedSupplementaryIndices())))
            {
                expanded.add(merged);
            }
        }

        return primaryAlignments.withLiftedAlignments(List.copyOf(expanded));
    }

    public Resolution finish(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords,
            final LiftedRecord primaryResult)
    {
        Set<Integer> absorbed = absorbedSupplementaries(primaryResult);
        if(hasSupplementaryMerge(primaryResult))
        {
            ++mStats.SupplementaryMerges;
            mStats.SupplementariesAbsorbed += absorbed.size();
        }

        List<LiftedRecord> selected = selectRecordAlignments(
                records, liftedRecords, primaryResult, absorbed);
        return new Resolution(annotateSpliceStrands(finaliseSupplementaries(records, selected, absorbed)), absorbed);
    }

    private static List<LiftedRecord> finaliseSupplementaries(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords,
            final Set<Integer> absorbedSupplementaries)
    {
        List<LiftedRecord> finalised = null;
        for(int i = 1; i < records.size(); ++i)
        {
            if(absorbedSupplementaries.contains(i))
            {
                continue;
            }

            LiftedRecord lifted = liftedRecords.get(i);
            LiftedRecord result = AlignmentSelector.finaliseSupplementary(lifted);
            if(result != lifted)
            {
                if(finalised == null)
                {
                    finalised = new ArrayList<>(liftedRecords);
                }
                finalised.set(i, result);
            }
        }
        return finalised != null ? List.copyOf(finalised) : liftedRecords;
    }

    private List<LiftedRecord> selectRecordAlignments(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords,
            final LiftedRecord primaryResult, final Set<Integer> absorbedSupplementaries)
    {
        if(!primaryResult.hasPlacement())
        {
            return liftedRecords;
        }

        List<AlignmentSelector.RecordAlignment> placements = placements(records, liftedRecords);
        if(placements.isEmpty())
        {
            return liftedRecords;
        }

        LiftedAlignment primary = primaryResult.primaryAlignment();
        List<AlignmentSelector.RecordAlignment> selected = AlignmentSelector.selectSupplementaryAlignments(
                primary.liftedPosition().Chromosome, primary.forwardStrand(), primary.liftedPosition().Position, primary.liftedCigar(),
                records.get(0).getReadBases(), placements);

        List<LiftedRecord> revised = null;
        for(AlignmentSelector.RecordAlignment supplementary : selected)
        {
            int recordIndex = supplementary.recordIndex();
            if(absorbedSupplementaries.contains(recordIndex))
            {
                continue;
            }

            int alignmentIndex = supplementary.alignmentIndex();
            LiftedRecord lifted = liftedRecords.get(recordIndex);
            if(alignmentIndex == lifted.primaryIndex())
            {
                continue;
            }

            LiftedAlignment selectedAlignment = lifted.liftedAlignments().get(alignmentIndex);
            int score = mAlignmentScorer.scoreRecord(selectedAlignment, records.get(recordIndex));
            LiftedRecord selectedRecord = lifted.withPrimaryIndex(alignmentIndex)
                    .withPrimaryAlignment(selectedAlignment.withGenomicScore(score));

            if(revised == null)
            {
                revised = new ArrayList<>(liftedRecords);
            }
            revised.set(recordIndex, selectedRecord);
        }

        return revised != null ? List.copyOf(revised) : liftedRecords;
    }

    private static List<AlignmentSelector.RecordAlignment> placements(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords)
    {
        List<AlignmentSelector.RecordAlignment> alignments = new ArrayList<>();
        for(int i = 1; i < records.size(); ++i)
        {
            SAMRecord record = records.get(i);
            LiftedRecord lifted = liftedRecords.get(i);
            if(!lifted.hasPlacement())
            {
                continue;
            }

            for(int alignmentIndex = 0; alignmentIndex < lifted.liftedAlignments().size(); ++alignmentIndex)
            {
                LiftedAlignment alignment = lifted.liftedAlignments().get(alignmentIndex);
                if(alignment.dropped())
                {
                    continue;
                }

                alignments.add(new AlignmentSelector.RecordAlignment(
                        i, alignmentIndex, alignment, record.getMappingQuality()));
            }
        }
        return List.copyOf(alignments);
    }

    private List<LiftedRecord> annotateSpliceStrands(final List<LiftedRecord> liftedRecords)
    {
        List<LiftedRecord> annotated = null;
        for(int i = 0; i < liftedRecords.size(); ++i)
        {
            LiftedRecord lifted = liftedRecords.get(i);
            if(!lifted.hasPlacement() || !lifted.hasNCigar() || lifted.transcriptStrand() != null)
            {
                continue;
            }

            Orientation strand = mMerger.spliceStrand(
                    lifted.finalChromosome(), lifted.finalPos(), lifted.finalCigar());
            if(strand != null)
            {
                if(annotated == null)
                {
                    annotated = new ArrayList<>(liftedRecords);
                }
                annotated.set(i, lifted.withPrimaryAlignment(lifted.primaryAlignment().withTranscriptStrand(strand)));
            }
        }
        return annotated != null ? List.copyOf(annotated) : liftedRecords;
    }

    private static Set<Integer> absorbedSupplementaries(final LiftedRecord primaryResult)
    {
        if(!hasSupplementaryMerge(primaryResult))
        {
            return Set.of();
        }

        return Set.copyOf(primaryResult.primaryAlignment().mergedSupplementaryIndices());
    }

    private static boolean hasSupplementaryMerge(final LiftedRecord primaryResult)
    {
        return primaryResult.hasPlacement()
                && primaryResult.primaryAlignment().hasSupplementaryMerge();
    }

    public record Resolution(
            List<LiftedRecord> liftedRecords, Set<Integer> absorbedSupplementaries)
    {
    }

    private record SupportedAlignmentKey(
            AlignmentKey alignment, List<Integer> absorbedSupplementaries)
    {
        private SupportedAlignmentKey
        {
            absorbedSupplementaries = List.copyOf(absorbedSupplementaries);
        }
    }

}
