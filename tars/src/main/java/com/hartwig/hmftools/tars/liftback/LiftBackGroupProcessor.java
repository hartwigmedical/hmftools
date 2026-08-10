package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;
import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.features.GenomicAlignmentScorer;
import com.hartwig.hmftools.tars.liftback.features.OverhangGate;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver;

import htsjdk.samtools.SAMRecord;

public class LiftBackGroupProcessor
{
    private final LiftBackDiscriminator mDiscriminator;
    private final SupplementaryResolver mSupplementaryResolver;
    private final OverhangGate mOverhangGate;
    private final GenomicAlignmentScorer mAlignmentScorer;
    private final BamRecordEmitter mEmitter;

    private int mRecordsSeen;
    private int mPrimariesUnmappedExcludedRegion;
    private int mPrimariesUnmappedOverCap;
    private int mPrimariesUnmappedLowAlignmentScore;
    private int mSupplementaryCandidates;
    private int mPrimaryRevisions;
    private int mSupplementaryMerges;
    private int mSupplementariesAbsorbed;

    public LiftBackGroupProcessor(
            final LiftBackDiscriminator discriminator, final SupplementaryResolver supplementaryResolver,
            final OverhangGate overhangGate, final GenomicAlignmentScorer alignmentScorer,
            final RefGenomeInterface refGenome, final ExcludedRegions excludedRegions)
    {
        mDiscriminator = discriminator;
        mSupplementaryResolver = supplementaryResolver;
        mOverhangGate = overhangGate;
        mAlignmentScorer = alignmentScorer;
        mEmitter = new BamRecordEmitter(
                discriminator.contigTranslator(), supplementaryResolver != null, refGenome, excludedRegions);
    }

    public int recordsSeen() { return mRecordsSeen; }

    public int primariesSeen() { return mEmitter.primariesSeen(); }

    public int primariesLiftFailed() { return mEmitter.primariesLiftFailed(); }

    public int primariesUnmappedExcludedRegion() { return mPrimariesUnmappedExcludedRegion; }

    public int primariesUnmappedOverCap() { return mPrimariesUnmappedOverCap; }

    public int primariesUnmappedLowAlignmentScore() { return mPrimariesUnmappedLowAlignmentScore; }

    public int supplementaryCandidates() { return mSupplementaryCandidates; }

    public int primaryRevisions() { return mPrimaryRevisions; }

    public int supplementaryMerges() { return mSupplementaryMerges; }

    public int supplementariesAbsorbed() { return mSupplementariesAbsorbed; }

    public void processNameGroup(final List<SAMRecord> group, final Consumer<SAMRecord> consumer)
    {
        mRecordsSeen += group.size();

        List<SAMRecord> firstOfPair = new ArrayList<>();
        List<SAMRecord> secondOfPair = new ArrayList<>();
        for(SAMRecord record : group)
        {
            (firstInPair(record) ? firstOfPair : secondOfPair).add(record);
        }

        LiftedRecord firstAlignments = liftPrimaryAlignments(firstOfPair);
        LiftedRecord secondAlignments = liftPrimaryAlignments(secondOfPair);
        LiftedMatePair matePair = new LiftedMatePair();

        ReadDecision firstDecision = decideRead(
                firstOfPair, List.of(), secondAlignments, firstAlignments);
        recordPrimary(firstOfPair, firstDecision, matePair);

        ReadDecision secondDecision = decideRead(
                secondOfPair, firstDecision.introducedIntrons(),
                chosenPlacement(firstDecision.primaryResult(), firstAlignments), secondAlignments);
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
        return mDiscriminator.liftPrimaryAlignments(records.get(0), mOverhangGate);
    }

    private ReadDecision decideRead(
            final List<SAMRecord> records, final List<ChrBaseRegion> mateHintIntrons,
            final LiftedRecord mate, final LiftedRecord preLiftedPrimaryAlignments)
    {
        if(records.isEmpty())
        {
            return ReadDecision.EMPTY;
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
                ? preLiftedPrimaryAlignments : mDiscriminator.liftPrimaryAlignments(primary, mOverhangGate);
        boolean hasSupplementaries = records.size() > 1;

        if(hasSupplementaries)
        {
            // Step 2: lift the supplementary records and add any supported spliced primary alignments.
            for(int i = 1; i < records.size(); ++i)
            {
                resolved.set(i, mDiscriminator.liftSupplementaryAlignment(records.get(i)));
            }
            primaryAlignments = addSupplementarySupportedAlignments(
                    records, primaryAlignments, resolved, mateHintIntrons);
        }

        // Step 3: score the genomic alignments and select the primary placement.
        resolved.set(0, selectPrimaryAlignment(primary, primaryAlignments, mate));
        SupplementaryResolution supplementary = supplementaryResult(records.size(), resolved.get(0));
        List<ChrBaseRegion> introducedIntrons = new ArrayList<>();

        if(supplementary.revised())
        {
            LiftedAlignment winner = resolved.get(0).primaryAlignment();
            introducedIntrons.addAll(winner.MergedSupplementaryIntrons);
            ++mPrimaryRevisions;
            ++mSupplementaryMerges;
            mSupplementariesAbsorbed += supplementary.absorbed().size();
        }

        List<LiftedRecord> finalRecords = annotateSpliceStrands(List.copyOf(resolved));

        UnmapDecision unmapDecision = applyUnmapPolicy(primary, finalRecords.get(0));
        if(unmapDecision.result() != finalRecords.get(0))
        {
            finalRecords = replacePrimary(finalRecords, unmapDecision.result());
        }
        return new ReadDecision(
                finalRecords, supplementary.absorbed(), introducedIntrons, unmapDecision.unmapped());
    }

    private LiftedRecord addSupplementarySupportedAlignments(
            final List<SAMRecord> records, final LiftedRecord primaryAlignments,
            final List<LiftedRecord> liftedRecords, final List<ChrBaseRegion> mateHintIntrons)
    {
        SAMRecord primary = records.get(0);
        if(mSupplementaryResolver == null || primary.getReadUnmappedFlag() || !primaryAlignments.hasPlacement())
        {
            return primaryAlignments;
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
            return primaryAlignments;
        }
        ++mSupplementaryCandidates;

        List<LiftedAlignment> expanded = new ArrayList<>(primaryAlignments.liftedAlignments());
        Set<SupplementaryAlignmentKey> seen = new HashSet<>();
        for(LiftedAlignment alignment : primaryAlignments.liftedAlignments())
        {
            seen.add(new SupplementaryAlignmentKey(alignment.key(), List.of()));
        }

        for(LiftedAlignment alignment : primaryAlignments.liftedAlignments())
        {
            if(alignment.Dropped || alignment.LiftedCigar == null)
            {
                continue;
            }

            SupplementaryResolver.Result result = mSupplementaryResolver.resolve(new SupplementaryResolver.Candidate(
                    alignment.LiftedChromosome, alignment.ForwardStrand, primary.getReadLength(),
                    alignment.LiftedPos, alignment.LiftedCigar, supplementaries,
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

            LiftedAlignment merged = alignment.withSupplementaryMerge(
                    result.mergedStart(), result.mergedCigar(), result.droppedSupplementaryIndices(),
                    result.introducedIntrons(), supplementaryMapQuality, result.spliceStrand());
            if(seen.add(new SupplementaryAlignmentKey(merged.key(), result.droppedSupplementaryIndices())))
            {
                expanded.add(merged);
            }
        }

        return primaryAlignments.withLiftedAlignments(List.copyOf(expanded));
    }

    private LiftedRecord selectPrimaryAlignment(
            final SAMRecord primary, final LiftedRecord alignments, final LiftedRecord mate)
    {
        if(!alignments.hasPlacement())
        {
            return alignments;
        }
        if(mAlignmentScorer != null)
        {
            mAlignmentScorer.scoreCandidates(alignments.liftedAlignments(), primary);
        }
        return mDiscriminator.selectPrimaryAlignment(primary, alignments.liftedAlignments(), mate);
    }

    private List<LiftedRecord> annotateSpliceStrands(final List<LiftedRecord> liftedRecords)
    {
        if(mSupplementaryResolver == null)
        {
            return liftedRecords;
        }

        List<LiftedRecord> annotated = new ArrayList<>(liftedRecords);
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
                annotated.set(i, lifted.withPrimaryTranscriptStrand(strand));
            }
        }
        return List.copyOf(annotated);
    }

    private Integer finalAlignmentScore(final SAMRecord primary, final LiftedRecord lifted)
    {
        // Lift-only test paths omit supplementary resolution and intentionally skip the production AS floor.
        if(mSupplementaryResolver == null)
        {
            return null;
        }
        if(lifted.hasPlacement() && mAlignmentScorer != null)
        {
            int score = mAlignmentScorer.scoreRecord(lifted.primaryAlignment(), primary);
            if(score != Integer.MIN_VALUE)
            {
                return score;
            }
        }
        return primary.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
    }

    private UnmapDecision applyUnmapPolicy(final SAMRecord primary, final LiftedRecord result)
    {
        if(mEmitter.excludes(result))
        {
            ++mPrimariesUnmappedExcludedRegion;
            return new UnmapDecision(LiftedRecord.unmapped("excluded_region_unmapped"), true);
        }
        if(primary.getReadUnmappedFlag())
        {
            return new UnmapDecision(result, false);
        }
        if(BamRecordEmitter.exceedsMappingCap(primary, result))
        {
            ++mPrimariesUnmappedOverCap;
            TARS_LOGGER.trace("over-cap unmap {}:{}: inputMapQuality=0, no XA",
                    primary.getReferenceName(), primary.getAlignmentStart());
            return new UnmapDecision(LiftedRecord.unmapped("over_cap_unmapped"), true);
        }

        Integer alignmentScore = finalAlignmentScore(primary, result);
        if(alignmentScore != null && alignmentScore < PRIMARY_AS_UNMAP_THRESHOLD)
        {
            ++mPrimariesUnmappedLowAlignmentScore;
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

    private static SupplementaryResolution supplementaryResult(
            final int recordCount, final LiftedRecord primaryResult)
    {
        if(primaryResult == null || !primaryResult.hasPlacement()
                || !primaryResult.primaryAlignment().hasSupplementaryMerge())
        {
            return SupplementaryResolution.NONE;
        }

        Set<Integer> absorbed = new HashSet<>();
        for(Integer index : primaryResult.primaryAlignment().MergedSupplementaryIndices)
        {
            if(index >= 0 && index < recordCount)
            {
                absorbed.add(index);
            }
        }
        return new SupplementaryResolution(absorbed, true);
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

    private record SupplementaryResolution(Set<Integer> absorbed, boolean revised)
    {
        private static final SupplementaryResolution NONE = new SupplementaryResolution(Set.of(), false);

        private SupplementaryResolution
        {
            absorbed = Set.copyOf(absorbed);
        }
    }

    private record UnmapDecision(LiftedRecord result, boolean unmapped)
    {
    }

    private record SupplementaryAlignmentKey(AlignmentKey alignment, List<Integer> absorbedSupplementaries)
    {
        private SupplementaryAlignmentKey
        {
            absorbedSupplementaries = List.copyOf(absorbedSupplementaries);
        }
    }
}
