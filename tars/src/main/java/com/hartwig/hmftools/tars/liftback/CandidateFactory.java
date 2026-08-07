package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.tars.common.TarsConstants.CONFIDENT_MAPQ;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.features.OverhangGate;
import com.hartwig.hmftools.tars.liftback.features.SoftClipExtender;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver;

import htsjdk.samtools.SAMRecord;

public final class CandidateFactory
{
    public record CandidateSet(LiftedRecord record)
    {
        public List<LiftedAlignment> alignments()
        {
            return record.liftedAlignments();
        }

        public boolean hasPlacement()
        {
            return record.hasPlacement();
        }
    }

    public record Expansion(CandidateSet candidates, boolean attempted)
    {
    }

    public record Revision(LiftedRecord result, boolean revised, List<ChrBaseRegion> introducedIntrons)
    {
        public Revision
        {
            introducedIntrons = List.copyOf(introducedIntrons);
        }
    }

    private record MergeCandidateKey(AlignmentKey alignment, List<Integer> absorbedSupplementaries)
    {
        private MergeCandidateKey
        {
            absorbedSupplementaries = List.copyOf(absorbedSupplementaries);
        }
    }

    private final LiftBackDiscriminator mDiscriminator;
    private final SupplementaryResolver mSupplementaryResolver;
    private final OverhangGate mOverhangGate;
    private final SoftClipExtender mSoftClipExtender;

    public CandidateFactory(
            final LiftBackDiscriminator discriminator, final SupplementaryResolver supplementaryResolver,
            final OverhangGate overhangGate, final SoftClipExtender softClipExtender)
    {
        mDiscriminator = discriminator;
        mSupplementaryResolver = supplementaryResolver;
        mOverhangGate = overhangGate;
        mSoftClipExtender = softClipExtender;
    }

    public CandidateSet liftPrimary(final SAMRecord primary)
    {
        return new CandidateSet(mDiscriminator.liftPrimaryCandidates(primary, mOverhangGate));
    }

    public LiftedRecord liftSupplementary(final SAMRecord supplementary)
    {
        return mDiscriminator.resolve(supplementary, mOverhangGate, null, null);
    }

    public LiftedRecord selectPrimary(
            final SAMRecord primary, final CandidateSet candidates, final LiftedRecord mate)
    {
        if(!candidates.hasPlacement())
        {
            return candidates.record();
        }
        return mDiscriminator.resolvePrimaryCandidates(
                primary, candidates.alignments(), mSoftClipExtender, mate);
    }

    public Expansion addSupplementaryCandidates(
            final List<SAMRecord> records, final CandidateSet primaryCandidates,
            final List<LiftedRecord> liftedRecords, final List<ChrBaseRegion> mateHintIntrons)
    {
        SAMRecord primary = records.get(0);
        if(mSupplementaryResolver == null || primary.getReadUnmappedFlag() || !primaryCandidates.hasPlacement())
        {
            return new Expansion(primaryCandidates, false);
        }

        List<SupplementaryResolver.Supplementary> supplementaries = new ArrayList<>();
        for(int i = 1; i < records.size(); ++i)
        {
            SAMRecord record = records.get(i);
            LiftedRecord lifted = liftedRecords.get(i);
            if(record.getSupplementaryAlignmentFlag() && !record.getReadUnmappedFlag()
                    && lifted != null && lifted.hasPlacement())
            {
                supplementaries.add(new SupplementaryResolver.Supplementary(
                        i, lifted.finalChromosome(), !record.getReadNegativeStrandFlag(),
                        lifted.finalPos(), lifted.finalCigar(), record.getMappingQuality()));
            }
        }
        if(supplementaries.isEmpty())
        {
            return new Expansion(primaryCandidates, false);
        }

        List<LiftedAlignment> expanded = new ArrayList<>(primaryCandidates.alignments());
        Set<MergeCandidateKey> seen = new HashSet<>();
        for(LiftedAlignment alignment : primaryCandidates.alignments())
        {
            seen.add(new MergeCandidateKey(alignment.key(), List.of()));
        }

        for(LiftedAlignment candidate : primaryCandidates.alignments())
        {
            if(candidate.Dropped || candidate.LiftedCigar == null)
            {
                continue;
            }

            SupplementaryResolver.Result result = mSupplementaryResolver.resolve(new SupplementaryResolver.Candidate(
                    candidate.LiftedChromosome, candidate.ForwardStrand, primary.getReadLength(),
                    candidate.LiftedPos, candidate.LiftedCigar, supplementaries,
                    primary.getReadBases(), mateHintIntrons));
            if(!result.merged())
            {
                continue;
            }

            int supplementaryMapQuality = 0;
            for(Integer recordIndex : result.droppedSupplementaryIndices())
            {
                supplementaryMapQuality = Math.max(
                        supplementaryMapQuality, records.get(recordIndex).getMappingQuality());
            }

            LiftedAlignment merged = candidate.withSupplementaryMerge(
                    result.mergedStart(), result.mergedCigar(), result.droppedSupplementaryIndices(),
                    result.introducedIntrons(), supplementaryMapQuality, result.spliceStrand());
            if(seen.add(new MergeCandidateKey(merged.key(), result.droppedSupplementaryIndices())))
            {
                expanded.add(merged);
            }
        }

        LiftedRecord expandedRecord = primaryCandidates.record().withLiftedAlignments(List.copyOf(expanded));
        return new Expansion(new CandidateSet(expandedRecord), true);
    }

    public Revision reviseWithoutSupplementaries(
            final SAMRecord primary, final LiftedRecord lifted, final List<ChrBaseRegion> mateHintIntrons)
    {
        if(mSupplementaryResolver == null || primary.getReadUnmappedFlag() || !lifted.hasPlacement())
        {
            return new Revision(lifted, false, List.of());
        }

        SupplementaryResolver.Result result = mSupplementaryResolver.resolve(new SupplementaryResolver.Candidate(
                lifted.finalChromosome(), !primary.getReadNegativeStrandFlag(), primary.getReadLength(),
                lifted.finalPos(), lifted.finalCigar(), List.of(), primary.getReadBases(), mateHintIntrons));
        if(!result.merged())
        {
            return new Revision(lifted, false, List.of());
        }

        int mapQuality = lifted.numLoci() == 1 ? CONFIDENT_MAPQ : lifted.updatedMapQuality();
        LiftedRecord revised = lifted.withRevisedPrimary(
                result.mergedStart(), result.mergedCigar(), mapQuality, "supplementary-resolved");
        TARS_LOGGER.trace(
                "supplementary-resolved {}: {}:{} {} -> {} (depth {})",
                primary.getReadName(), lifted.finalChromosome(), lifted.finalPos(), lifted.finalCigar(),
                result.mergedCigar(), result.chainDepth());
        return new Revision(revised, true, result.introducedIntrons());
    }

    public List<LiftedRecord> extendSoftClips(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords)
    {
        if(mSoftClipExtender == null)
        {
            return liftedRecords;
        }
        List<LiftedRecord> extended = new ArrayList<>(liftedRecords);
        for(int i = 0; i < liftedRecords.size(); ++i)
        {
            LiftedRecord lifted = liftedRecords.get(i);
            if(lifted != null && lifted.hasPlacement())
            {
                extended.set(i, lifted.withExtendedSoftClips(mSoftClipExtender, records.get(i)));
            }
        }
        return List.copyOf(extended);
    }

    public List<LiftedRecord> snapBoundaries(
            final List<SAMRecord> records, final List<LiftedRecord> liftedRecords)
    {
        if(mSupplementaryResolver == null || records.size() < 2)
        {
            return liftedRecords;
        }
        List<LiftedRecord> snapped = new ArrayList<>(liftedRecords);
        for(int i = 0; i < records.size(); ++i)
        {
            LiftedRecord lifted = liftedRecords.get(i);
            if(lifted == null || !lifted.hasPlacement())
            {
                continue;
            }
            SupplementaryResolver.BoundarySnap snap = mSupplementaryResolver.snapToAnnotatedBoundary(
                    lifted.finalChromosome(), lifted.finalPos(), lifted.finalCigar());
            if(snap != null)
            {
                snapped.set(i, lifted.withRevisedPrimary(
                        snap.start(), snap.cigar(), lifted.updatedMapQuality(), "boundary-snapped"));
            }
        }
        return List.copyOf(snapped);
    }

    public List<LiftedRecord> assignSpliceStrands(final List<LiftedRecord> liftedRecords)
    {
        if(mSupplementaryResolver == null)
        {
            return liftedRecords;
        }
        List<LiftedRecord> assigned = new ArrayList<>(liftedRecords);
        for(int i = 0; i < liftedRecords.size(); ++i)
        {
            LiftedRecord lifted = liftedRecords.get(i);
            if(lifted == null || !lifted.hasPlacement() || !lifted.hasNCigar() || lifted.transcriptStrand() != 0)
            {
                continue;
            }
            int strand = mSupplementaryResolver.spliceStrand(
                    lifted.finalChromosome(), lifted.finalPos(), lifted.finalCigar());
            if(strand != 0)
            {
                assigned.set(i, lifted.withPrimaryTranscriptStrand(strand));
            }
        }
        return List.copyOf(assigned);
    }

    public Integer finalAlignmentScore(final SAMRecord primary, final LiftedRecord lifted)
    {
        if(mSupplementaryResolver == null)
        {
            return null;
        }
        if(lifted.hasPlacement() && mSoftClipExtender != null)
        {
            int score = mSoftClipExtender.scoreRecord(lifted.primaryAlignment(), primary);
            if(score != Integer.MIN_VALUE)
            {
                return score;
            }
        }
        return primary.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
    }
}
