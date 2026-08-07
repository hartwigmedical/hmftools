package com.hartwig.hmftools.tars.liftback;

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
import com.hartwig.hmftools.tars.liftback.CandidateFactory.CandidateSet;
import com.hartwig.hmftools.tars.liftback.features.OverhangGate;
import com.hartwig.hmftools.tars.liftback.features.SoftClipExtender;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver;

import htsjdk.samtools.SAMRecord;

public class LiftBackGroupProcessor
{
    private final CandidateFactory mCandidateFactory;
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
            final OverhangGate overhangGate, final SoftClipExtender softClipExtender,
            final RefGenomeInterface refGenome, final ExcludedRegions excludedRegions)
    {
        mCandidateFactory = new CandidateFactory(
                discriminator, supplementaryResolver, overhangGate, softClipExtender);
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

        CandidateSet firstCandidates = liftCandidates(firstOfPair);
        CandidateSet secondCandidates = liftCandidates(secondOfPair);
        LiftedMatePair matePair = new LiftedMatePair();

        ReadDecision firstDecision = decideRead(
                firstOfPair, List.of(), candidateRecord(secondCandidates), firstCandidates);
        recordPrimary(firstOfPair, firstDecision, matePair);

        ReadDecision secondDecision = decideRead(
                secondOfPair, firstDecision.introducedIntrons(),
                chosenPlacement(firstDecision.primaryResult(), candidateRecord(firstCandidates)), secondCandidates);
        recordPrimary(secondOfPair, secondDecision, matePair);

        emit(firstOfPair, firstDecision, matePair, consumer);
        emit(secondOfPair, secondDecision, matePair, consumer);
    }

    private CandidateSet liftCandidates(final List<SAMRecord> records)
    {
        if(records.isEmpty() || records.get(0).getReadUnmappedFlag())
        {
            return null;
        }
        return mCandidateFactory.liftPrimary(records.get(0));
    }

    private ReadDecision decideRead(
            final List<SAMRecord> records, final List<ChrBaseRegion> mateHintIntrons,
            final LiftedRecord mate, final CandidateSet liftedPrimaryCandidates)
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
        CandidateSet primaryCandidates = liftedPrimaryCandidates != null
                ? liftedPrimaryCandidates : mCandidateFactory.liftPrimary(primary);
        boolean hasSupplementaries = records.size() > 1;

        if(hasSupplementaries)
        {
            for(int i = 1; i < records.size(); ++i)
            {
                resolved.set(i, mCandidateFactory.liftSupplementary(records.get(i)));
            }
            CandidateFactory.Expansion expansion = mCandidateFactory.addSupplementaryCandidates(
                    records, primaryCandidates, resolved, mateHintIntrons);
            if(expansion.attempted())
            {
                ++mSupplementaryCandidates;
            }
            primaryCandidates = expansion.candidates();
        }

        resolved.set(0, mCandidateFactory.selectPrimary(primary, primaryCandidates, mate));
        SupplementaryResolution supplementary = supplementaryResult(records.size(), resolved.get(0));
        List<ChrBaseRegion> introducedIntrons = new ArrayList<>();

        if(!hasSupplementaries)
        {
            CandidateFactory.Revision revision = mCandidateFactory.reviseWithoutSupplementaries(
                    primary, resolved.get(0), mateHintIntrons);
            resolved.set(0, revision.result());
            if(revision.revised())
            {
                supplementary = new SupplementaryResolution(Set.of(), true);
                introducedIntrons.addAll(revision.introducedIntrons());
                ++mPrimaryRevisions;
            }
        }
        else if(supplementary.revised())
        {
            LiftedAlignment winner = resolved.get(0).primaryAlignment();
            introducedIntrons.addAll(winner.MergedSupplementaryIntrons);
            ++mPrimaryRevisions;
            ++mSupplementaryMerges;
            mSupplementariesAbsorbed += supplementary.absorbed().size();
        }

        List<LiftedRecord> finalRecords = List.copyOf(resolved);
        if(hasSupplementaries)
        {
            finalRecords = mCandidateFactory.extendSoftClips(records, finalRecords);
        }
        if(!supplementary.revised())
        {
            finalRecords = mCandidateFactory.snapBoundaries(records, finalRecords);
        }
        finalRecords = mCandidateFactory.assignSpliceStrands(finalRecords);

        UnmapDecision unmapDecision = applyUnmapPolicy(primary, finalRecords.get(0));
        if(unmapDecision.result() != finalRecords.get(0))
        {
            finalRecords = replacePrimary(finalRecords, unmapDecision.result());
        }
        return new ReadDecision(
                finalRecords, supplementary.absorbed(), introducedIntrons, unmapDecision.unmapped());
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
            TARS_LOGGER.trace("over-cap unmap {}: inputMapQuality=0, no XA", primary.getReadName());
            return new UnmapDecision(LiftedRecord.unmapped("over_cap_unmapped"), true);
        }

        Integer alignmentScore = mCandidateFactory.finalAlignmentScore(primary, result);
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

    private static LiftedRecord candidateRecord(final CandidateSet candidates)
    {
        return candidates != null ? candidates.record() : null;
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
}
