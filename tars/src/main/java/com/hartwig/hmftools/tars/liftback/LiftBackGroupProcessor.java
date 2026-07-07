package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.LiftBackRecordOps.markPrimaryUnmapped;
import static com.hartwig.hmftools.tars.liftback.LiftBackRecordOps.refreshNmDropMd;
import static com.hartwig.hmftools.tars.liftback.MateFieldPatcher.patchMateFields;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.rewriteSaTag;
import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.CONFIDENT_MAPQ;
import static com.hartwig.hmftools.tars.common.TarsConstants.MAX_CONFIDENT_NM_FRACTION;
import static com.hartwig.hmftools.tars.common.TarsConstants.SUPP_AS_DROP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGate;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGateStatistics;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryResolver;
import com.hartwig.hmftools.tars.liftback.supplementary.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryCandidate;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryResult;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryStatistics;
import com.hartwig.hmftools.tars.liftback.supplementary.SupplementaryRecord;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;

// Stateless-per-group transform engine: resolves one contiguous read-name group to genomic coordinates,
// runs the optional supplementary-resolve / overhang-gate passes, patches mate fields against a
// provided LiftedMateInfoCache, and emits each kept record to an EmitSink. The cache is supplied by the
// caller so the standalone TarsApplication can pass its whole-sample pass-1 cache while the REDUX worker
// passes a fresh per-group cache (a complete name-sorted group is self-sufficient -- see
// [[project_redux_migration_scoping]]). Extracted verbatim from TarsApplication so both paths share one
// implementation.
public class LiftBackGroupProcessor
{
    private final LiftBackResolver mResolver;

    // optional supplementary resolver: merges primary + supp across annotated junctions. Null when disabled.
    private final SupplementaryResolver mSupplementaryResolver;

    // null when no ref genome is loaded. Runs the overhang gate: iterative terminal collapse + softclip reclaim.
    private final OverhangGate mOverhangGate;

    // genomic reference for the post-lift NM recompute. Null when no ref is loaded (supplementary resolve + extend
    // both off), in which case NM is cleared rather than recomputed.
    private final RefSequenceSource mRefSource;

    // genomic excluded regions (rRNA / contamination), checked post-lift. Null when no regions file configured.
    private final ExcludedRegions mExcludedRegions;

    private final LiftBackStats mStats;

    // primaries unmapped because bwa exceeded the XA cap (MAPQ 0 + no XA = maps to too many loci).
    private long mOverCapUnmapped;

    // reads dropped because they (primary -> whole fragment, or a supp) lift into an excluded region.
    private long mExcludedReads;

    public long overCapUnmapped() { return mOverCapUnmapped; }

    public long excludedReads() { return mExcludedReads; }

    private static final String AS_TAG = "AS";

    // sink for emitted records: the standalone writes a BAM + TSV, the REDUX worker writes the shared BAM.
    public interface EmitSink
    {
        void emit(SAMRecord record, LiftBackResult result);
    }

    public LiftBackGroupProcessor(
            final LiftBackResolver resolver, final SupplementaryResolver supplementaryResolver,
            final OverhangGate overhangGate, final RefSequenceSource refSource,
            final ExcludedRegions excludedRegions, final LiftBackStats stats)
    {
        mResolver = resolver;
        mSupplementaryResolver = supplementaryResolver;
        mOverhangGate = overhangGate;
        mRefSource = refSource;
        mExcludedRegions = excludedRegions;
        mStats = stats;
    }

    // Processes all records sharing one read name as a group so primary + supps resolve together. This
    // lets supplementary resolve see all split-read components and lets /2 hint /1's introns (and vice versa).
    public void processNameGroup(
            final List<SAMRecord> group, final LiftedMateInfoCache liftedMateInfoCache, final EmitSink sink)
    {
        List<SAMRecord> firstOfPair = new ArrayList<>();
        List<SAMRecord> secondOfPair = new ArrayList<>();
        for(final SAMRecord record : group)
        {
            if(!record.getReadPairedFlag() || record.getFirstOfPairFlag())
            {
                firstOfPair.add(record);
            }
            else
            {
                secondOfPair.add(record);
            }
        }

        // Two passes per pair so mate /2 can use mate /1's resolved junctions as a hint, and a
        // re-decide of mate /1 picks up mate /2's hints when mate /1 wasn't initially resolved.
        PassCounterSnapshot beforeM1 = snapshotPassCounters();
        MateDecision m1 = decideMateGroup(firstOfPair, Collections.emptyList());
        PassCounterSnapshot afterM1 = snapshotPassCounters();
        refreshMateInfoCache(firstOfPair, m1, liftedMateInfoCache);
        MateDecision m2 = decideMateGroup(secondOfPair, m1.IntroducedIntrons);
        refreshMateInfoCache(secondOfPair, m2, liftedMateInfoCache);

        MateDecision m1Final;
        if(m1.IntroducedIntrons.isEmpty() && !m2.IntroducedIntrons.isEmpty())
        {
            // The provisional m1 decision is discarded and re-made with mate /2's introns. Roll back the engine
            // pass-effect counters it bumped (overhang-gate / supplementary resolve / excluded reads) so
            // the discarded pass is not double-counted; the re-decision below re-counts. Per-record stats are
            // unaffected (recorded once in writeMateGroup on the chosen decision); the emitted BAM is unaffected.
            rewindProvisionalCounters(beforeM1, afterM1);
            m1Final = decideMateGroup(firstOfPair, m2.IntroducedIntrons);
            refreshMateInfoCache(firstOfPair, m1Final, liftedMateInfoCache);
        }
        else
        {
            m1Final = m1;
        }

        writeMateGroup(firstOfPair, m1Final, liftedMateInfoCache, sink);
        writeMateGroup(secondOfPair, m2, liftedMateInfoCache, sink);
    }

    // Push the post-supplementary-resolve primary back into the mate cache so the partner mate's MC tag uses
    // the merged/extended cigar instead of the pre-resolve one the seed used.
    private void refreshMateInfoCache(
            final List<SAMRecord> mateRecords, final MateDecision decision, final LiftedMateInfoCache cache)
    {
        if(mateRecords.isEmpty() || decision.PrimaryResult == null)
        {
            return;
        }
        SAMRecord primary = null;
        for(final SAMRecord r : mateRecords)
        {
            if(!r.getSupplementaryAlignmentFlag())
            {
                primary = r;
                break;
            }
        }
        if(primary == null)
        {
            return;
        }
        LiftedMateInfo refreshed = LiftBackRecordOps.toLiftedMateInfo(primary, decision.PrimaryResult);
        boolean firstOfPair = !primary.getReadPairedFlag() || primary.getFirstOfPairFlag();
        cache.recordPrimaryAlignment(primary.getReadName(), firstOfPair, refreshed);
    }

    // Snapshot of the pass-effect counters mutated by one decideMateGroup pass, so a discarded provisional mate
    // decision can be rolled back rather than double-counted in the summary (see processNameGroup).
    private static final class PassCounterSnapshot
    {
        final long[] Overhang;
        final int[] SupplementaryResolve;
        final long ExcludedReads;

        PassCounterSnapshot(final long[] overhang, final int[] supplementaryResolve, final long excludedReads)
        {
            Overhang = overhang;
            SupplementaryResolve = supplementaryResolve;
            ExcludedReads = excludedReads;
        }
    }

    private PassCounterSnapshot snapshotPassCounters()
    {
        return new PassCounterSnapshot(
                mOverhangGate != null ? mOverhangGate.statistics().snapshot() : null,
                mSupplementaryResolver != null ? mSupplementaryResolver.statistics().snapshot() : null,
                mExcludedReads);
    }

    // Subtract the (after - before) delta a discarded provisional pass contributed, leaving only the kept counts.
    private void rewindProvisionalCounters(final PassCounterSnapshot before, final PassCounterSnapshot after)
    {
        if(mOverhangGate != null)
        {
            OverhangGateStatistics overhang = mOverhangGate.statistics();
            long[] current = overhang.snapshot();
            for(int i = 0; i < current.length; ++i)
            {
                current[i] -= after.Overhang[i] - before.Overhang[i];
            }
            overhang.restore(current);
        }
        if(mSupplementaryResolver != null)
        {
            SupplementaryStatistics supplementary = mSupplementaryResolver.statistics();
            int[] current = supplementary.snapshot();
            for(int i = 0; i < current.length; ++i)
            {
                current[i] -= after.SupplementaryResolve[i] - before.SupplementaryResolve[i];
            }
            supplementary.restore(current);
        }
        mExcludedReads -= after.ExcludedReads - before.ExcludedReads;
    }

    // Per-mate decision: resolves every record, runs the supplementary resolver, returns the lifted results +
    // per-record drop flags + introns introduced by supplementary resolve (for mate hinting). Does not emit.
    private static final class MateDecision
    {
        LiftBackResult[] Resolved;
        boolean[] DroppedBySupplementary;
        LiftBackResult PrimaryResult;
        List<ChrBaseRegion> IntroducedIntrons;
        int PrimaryIndex;
        boolean PrimaryPostProcessed;

        MateDecision(final LiftBackResult[] resolved, final boolean[] droppedBySupplementary,
                final LiftBackResult primaryResult, final List<ChrBaseRegion> introducedIntrons,
                final int primaryIndex, final boolean primaryPostProcessed)
        {
            Resolved = resolved;
            DroppedBySupplementary = droppedBySupplementary;
            PrimaryResult = primaryResult;
            IntroducedIntrons = introducedIntrons != null ? introducedIntrons : Collections.emptyList();
            PrimaryIndex = primaryIndex;
            PrimaryPostProcessed = primaryPostProcessed;
        }

        static MateDecision empty()
        {
            return new MateDecision(new LiftBackResult[0], null, null, Collections.emptyList(), -1, false);
        }
    }

    private MateDecision decideMateGroup(final List<SAMRecord> records, final List<ChrBaseRegion> mateHintIntrons)
    {
        if(records.isEmpty())
        {
            return MateDecision.empty();
        }

        SAMRecord primary = null;
        for(final SAMRecord record : records)
        {
            if(record.getSupplementaryAlignmentFlag())
                continue;
            if(primary == null)
            {
                primary = record;
            }
        }

        // a valid BAM always has a primary in the mate group, so resolve it directly. Candidate cigars are peeled
        // by the overhang gate before the discriminator runs - see reconcileAlignmentsToGenome.
        LiftBackResult primaryResult = mResolver.resolve(primary, this::reconcileAlignmentsToGenome);

        LiftBackResult[] resolved = new LiftBackResult[records.size()];
        for(int i = 0; i < records.size(); i++)
        {
            SAMRecord record = records.get(i);
            resolved[i] = (record == primary && primaryResult != null)
                    ? primaryResult : mResolver.resolve(record, this::reconcileAlignmentsToGenome);
        }

        int primaryIdx = primary != null ? indexOfPrimary(records, primary) : -1;

        List<ChrBaseRegion> introduced = new ArrayList<>();
        boolean[] droppedBySupplementary = applySupplementaryResolve(records, resolved, primary, mateHintIntrons, introduced);

        // Finalize the chosen primary: peel it, then reclaim its terminal softclip (tx-match only). Not a redundant
        // repeat of the per-candidate peel in reconcileAlignmentsToGenome - it does real work in two cases:
        //  1. supplementary resolve merged a primary+supp into a fresh M N M cigar that never existed at candidate
        //     time, so it can carry its own fabricated terminal micro-junction; and
        //  2. the standalone softclip reclaim is deliberately deferred to here (post-resolve) so supplementary
        //     resolve saw the original clip - this is the chosen primary's first and only reclaim.
        // A clean primary is left unchanged by each step.
        reconcileChosenPrimary(records, resolved, primary, droppedBySupplementary);

        // post-lift exclusion: if the final primary placement lands in an excluded (rRNA / contamination) region,
        // unmap it REDUX-style - flip the result to UNMAPPED so the mate is coordinated via the cache (willBeUnmapped)
        // and the record is unmapped at apply time. Lifted genomic coords are the only space tx-contig reads can be
        // tested against the genomic region list, so this is post-lift, not the old pre-lift fragment pre-filter.
        if(primaryIdx >= 0 && liftsIntoExcludedRegion(resolved[primaryIdx]))
        {
            resolved[primaryIdx] = LiftBackResult.unmapped(LiftBackResult.RecordRole.PRIMARY, "excluded_region_unmapped");
            ++mExcludedReads;
        }

        // supplementary resolve / the overhang gate each replace resolved[primaryIdx] with a new result object when
        // they improve the primary. A changed reference means liftback post-processed it, so its stale bwa
        // AS must not be used to AS-filter it (see PRIMARY_AS_UNMAP_THRESHOLD).
        boolean primaryPostProcessed = primaryIdx >= 0 && resolved[primaryIdx] != primaryResult;

        return new MateDecision(resolved, droppedBySupplementary,
                primaryIdx >= 0 ? resolved[primaryIdx] : primaryResult,
                introduced, primaryIdx, primaryPostProcessed);
    }

    private static int indexOfPrimary(final List<SAMRecord> records, final SAMRecord primary)
    {
        for(int i = 0; i < records.size(); ++i)
            if(records.get(i) == primary)
            {
                return i;
            }
        return -1;
    }

    private void writeMateGroup(
            final List<SAMRecord> records, final MateDecision decision,
            final LiftedMateInfoCache liftedMateInfoCache, final EmitSink sink)
    {
        if(records.isEmpty())
        {
            return;
        }

        LiftBackResult[] resolved = decision.Resolved;
        boolean[] droppedBySupplementary = decision.DroppedBySupplementary;
        SAMRecord primary = null;
        for(SAMRecord r : records)
            if(!r.getSupplementaryAlignmentFlag())
            {
                primary = r;
                break;
            }

        // A supplementary whose primary ends up unmapped (e.g. the primary lifted into an excluded region) is
        // orphaned: there is no emitted primary for its SA to reference, so REDUX FragmentCoords would read a
        // dangling primary locus off it. Drop the whole group's supps in that case.
        boolean primaryUnmapped = decision.PrimaryResult == null || LiftBackRecordOps.willBeUnmapped(decision.PrimaryResult);

        // Dedup supplementaries that lift to the same (chrom, pos, cigar, strand) -- bwa can emit the
        // same junction across multiple transcript contigs and they collapse after liftback.
        Set<String> emittedSuppKeys = new HashSet<>();
        boolean[] willEmit = new boolean[records.size()];
        for(int i = 0; i < records.size(); i++)
        {
            SAMRecord record = records.get(i);
            LiftBackResult result = resolved[i];
            boolean drop = droppedBySupplementary != null && droppedBySupplementary[i];
            if(!drop && record.getSupplementaryAlignmentFlag() && primaryUnmapped)
            {
                drop = true;
                mStats.recordOrphanSuppDropped();
            }
            if(!drop && record != primary && record.getSupplementaryAlignmentFlag())
            {
                String key = dedupKey(result, record);
                if(!emittedSuppKeys.add(key))
                {
                    drop = true;
                }
            }
            // Drop supps the supplementary-resolve pass left behind that exist only because bwa-mem2 was run with
            // -T 19 below its default of 30 (see SUPP_AS_DROP_THRESHOLD). Only applies when supplementary resolve
            // ran, because a configuration without it might want to retain these supps for other reasons.
            if(!drop && mSupplementaryResolver != null && record.getSupplementaryAlignmentFlag()
                    && !record.getReadUnmappedFlag())
            {
                Integer alignmentScore = record.getIntegerAttribute(AS_TAG);
                if(alignmentScore != null && alignmentScore < SUPP_AS_DROP_THRESHOLD)
                {
                    drop = true;
                    mStats.recordLowAsSuppDropped();
                }
            }
            // a supp lifting into an excluded region is dropped (a supp can't be unmapped); its SA entry is
            // removed from the primary below so the primary doesn't reference a supp that isn't emitted.
            if(!drop && record.getSupplementaryAlignmentFlag() && liftsIntoExcludedRegion(result))
            {
                drop = true;
                ++mExcludedReads;
            }
            willEmit[i] = !drop;
        }

        // SA entries of dropped supps removed from the primary's SA so it never references a missing supp.
        Set<String> droppedSuppSaKeys = new HashSet<>();
        for(int i = 0; i < records.size(); i++)
        {
            if(willEmit[i])
                continue;
            SAMRecord record = records.get(i);
            LiftBackResult result = resolved[i];
            if(record == primary || !record.getSupplementaryAlignmentFlag() || result == null
                    || result.finalChrom() == null || result.finalCigar() == null
                    || result.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
                continue;
            droppedSuppSaKeys.add(suppSaKey(result));
        }

        // NH = the number of distinct genomic loci the read lifts back to, taken from the resolver's
        // chrom:pos-keyed locus count. Counting emitted records would inflate NH by tx-contig
        // representation multiplicity (one junction repeated across many transcript contigs all lifting
        // to the same locus).
        int nh = decision.PrimaryResult != null ? Math.max(decision.PrimaryResult.numLoci(), 1) : 1;

        for(int i = 0; i < records.size(); i++)
        {
            SAMRecord record = records.get(i);
            LiftBackResult result = resolved[i];
            if(!willEmit[i])
            {
                mStats.record(record, result);
                continue;
            }
            boolean primaryPostProcessed = i == decision.PrimaryIndex && decision.PrimaryPostProcessed;
            applyAndWriteRecord(record, result, nh, primaryPostProcessed, droppedSuppSaKeys, liftedMateInfoCache, sink);
        }
    }

    // builds a SupplementaryCandidate from the post-lift primary + its supplementary records and applies the
    // merge if accepted. Returns a parallel array marking which records were absorbed by the merge and
    // should be dropped from emission. Returns null when supplementary resolve is disabled / no-op.
    // introducedIntronsOut is appended-to (caller-provided list) so the partner mate can use them as hints in
    // the second pass.
    private boolean[] applySupplementaryResolve(
            final List<SAMRecord> records, final LiftBackResult[] resolved, final SAMRecord primary,
            final List<ChrBaseRegion> mateHintIntrons, final List<ChrBaseRegion> introducedIntronsOut)
    {
        if(mSupplementaryResolver == null || primary == null || primary.getReadUnmappedFlag())
        {
            return null;
        }

        int primaryIdx = -1;
        List<Integer> suppIndices = new ArrayList<>();
        for(int i = 0; i < records.size(); i++)
        {
            SAMRecord r = records.get(i);
            if(r == primary)
            {
                primaryIdx = i;
            }
            else if(r.getSupplementaryAlignmentFlag() && !r.getReadUnmappedFlag())
            {
                suppIndices.add(i);
            }
        }
        if(primaryIdx < 0)
        {
            return null;
        }

        LiftBackResult primaryRes = resolved[primaryIdx];
        if(primaryRes == null || primaryRes.finalCigar() == null
                || primaryRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
            return null;

        List<SupplementaryRecord> suppDtos = new ArrayList<>(suppIndices.size());
        for(int i = 0; i < suppIndices.size(); i++)
        {
            int idx = suppIndices.get(i);
            SAMRecord supp = records.get(idx);
            LiftBackResult suppRes = resolved[idx];
            if(suppRes == null || suppRes.finalCigar() == null
                    || suppRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
                continue;
            suppDtos.add(new SupplementaryRecord(
                    i, suppRes.finalChrom(), !supp.getReadNegativeStrandFlag(),
                    suppRes.finalPos(), suppRes.finalCigar(), supp.getMappingQuality()));
        }

        SupplementaryCandidate candidate = new SupplementaryCandidate(
                primaryRes.finalChrom(), !primary.getReadNegativeStrandFlag(), primary.getReadLength(),
                primaryRes.finalPos(), primaryRes.finalCigar(), primary.getMappingQuality(), suppDtos,
                primary.getReadBases(), mateHintIntrons);

        SupplementaryResult result = mSupplementaryResolver.resolve(candidate);
        if(!result.merged())
        {
            return null;
        }

        // rewrite the primary's LiftBackResult with the merged (now-spliced) cigar + start; the start may
        // shift for a left-extend. Merged MAPQ is the max of the primary and every merged supp; a merged
        // alignment left at 0 is bumped to 60 (the merge is the placement evidence). mark merged supps for drop.
        int mergedMapq = primaryRes.updatedMapq();
        for(Integer dtoIdx : result.droppedSupplementaryIndices())
        {
            mergedMapq = Math.max(mergedMapq, suppDtos.get(dtoIdx).mapq());
        }
        if(mergedMapq == 0)
        {
            mergedMapq = CONFIDENT_MAPQ;
        }
        resolved[primaryIdx] = primaryRes.withRevisedCigar(
                result.mergedStart(), result.mergedCigar(), true,
                mergedMapq, "supplementary-resolved");
        TARS_LOGGER.debug2("supplementary-resolved {}: {}:{} {} -> {} (depth {})",
                primary.getReadName(), primaryRes.finalChrom(), primaryRes.finalPos(),
                primaryRes.finalCigar(), result.mergedCigar(), result.chainDepth());

        if(introducedIntronsOut != null)
        {
            introducedIntronsOut.addAll(result.introducedIntrons());
        }

        boolean[] dropped = new boolean[records.size()];
        for(Integer dtoIdx : result.droppedSupplementaryIndices())
        {
            dropped[suppIndices.get(dtoIdx)] = true;
        }
        return dropped;
    }

    // Peels each candidate's terminal junction overhangs against the genome before the discriminator, so its
    // RefSoftClipped / RefFullMatch / TxHasNCigar inputs are measured facts, not raw-bwa guesses. Read bases are
    // taken in each candidate's own orientation; self and same-strand alts use the record bases as-is,
    // opposite-strand alts the reverse complement. Only the peel runs here (no standalone reclaim): self's softclip
    // is left intact for supplementary resolve (a primary+supp merge needs the clip), and the chosen primary is
    // reclaimed post-resolve in reconcileChosenPrimary. An XA alt the peel collapses to a purely contiguous alignment is a
    // fabricated placement and is dropped from the tag + the locus count; self is never dropped, only collapsed.
    private void reconcileAlignmentsToGenome(final List<LiftedAlignment> alignments, final SAMRecord record)
    {
        if(mOverhangGate == null || !mOverhangGate.enabled() || alignments.isEmpty())
        {
            return;
        }

        byte[] forwardBases = record.getReadBases();
        if(forwardBases == null || forwardBases.length == 0)
        {
            return;
        }

        LiftedAlignment self = alignments.get(0);
        boolean recordForward = !record.getReadNegativeStrandFlag();

        // Fast path for the dominant single-candidate read: peel self directly and skip the per-placement dedup map.
        // Self uses the record's own orientation (forward bases) and is never dropped.
        if(alignments.size() == 1)
        {
            if(self.LiftedCigar != null)
            {
                OverhangGate.Result peeled = mOverhangGate.peel(
                        self.LiftedChrom, self.LiftedPos, self.LiftedCigar, forwardBases);
                if(peeled.pos() != self.LiftedPos || !peeled.cigar().equals(self.LiftedCigar))
                {
                    alignments.set(0, self.withLiftedCigar(peeled.pos(), peeled.cigar()));
                }
            }
            return;
        }

        byte[] reverseBases = null;

        // Many candidates lift to an identical genomic placement (packed isoform contigs of one gene all collapsing
        // to one locus). peel() is deterministic in (chrom, pos, cigar, strand), so each distinct placement is
        // peeled once and the result reused for every duplicate.
        Map<String, OverhangGate.Result> placementsByKey = new HashMap<>();

        for(int i = 0; i < alignments.size(); ++i)
        {
            LiftedAlignment alignment = alignments.get(i);
            if(alignment.LiftedCigar == null)
            {
                continue;
            }

            String key = alignment.LiftedChrom + ':' + alignment.LiftedPos + ':' + alignment.LiftedCigar
                    + (alignment.ForwardStrand ? '+' : '-');
            OverhangGate.Result peeled = placementsByKey.get(key);
            if(peeled == null)
            {
                byte[] bases;
                if(alignment.ForwardStrand == recordForward)
                {
                    bases = forwardBases;
                }
                else
                {
                    if(reverseBases == null)
                    {
                        reverseBases = Arrays.copyOf(forwardBases, forwardBases.length);
                        SequenceUtil.reverseComplement(reverseBases);
                    }
                    bases = reverseBases;
                }

                peeled = mOverhangGate.peel(alignment.LiftedChrom, alignment.LiftedPos, alignment.LiftedCigar, bases);
                placementsByKey.put(key, peeled);
            }

            // An XA alt whose junction(s) fully collapsed is a fabricated placement: drop it from XA + the locus
            // count. Self keeps its collapsed cigar (it stays the record's placement).
            if(alignment != self && peeled.dropped())
            {
                alignment.Dropped = true;
                mOverhangGate.statistics().countAltDropped();
                continue;
            }

            if(peeled.pos() != alignment.LiftedPos || !peeled.cigar().equals(alignment.LiftedCigar))
            {
                alignments.set(i, alignment.withLiftedCigar(peeled.pos(), peeled.cigar()));
            }
        }
    }

    // True when the chosen primary's placement was lifted from a transcript contig; gates the standalone softclip
    // reclaim so a genomic/ref over-clip is left as bwa placed it.
    private static boolean isTxMatchPrimary(final LiftBackResult primaryRes)
    {
        for(final LiftedAlignment alignment : primaryRes.liftedAlignments())
        {
            if(alignment.IsPrimaryChoice)
            {
                return alignment.fromTxContig();
            }
        }
        return false;
    }

    // Finalize the chosen primary: peel its terminal junction overhangs (catches one merge fabricated in a fresh
    // M N M cigar), reclaim a surviving terminal softclip when the primary is a tx-match (deferred here so
    // supplementary resolve saw the original clip). Runs on whatever primary the discriminator or supplementary
    // resolve chose. A clean primary is unchanged by each step.
    private void reconcileChosenPrimary(
            final List<SAMRecord> records, final LiftBackResult[] resolved, final SAMRecord primary,
            final boolean[] droppedBySupplementary)
    {
        if(mOverhangGate == null)
        {
            return;
        }
        if(primary == null || primary.getReadUnmappedFlag())
        {
            return;
        }

        int primaryIdx = indexOfPrimary(records, primary);
        if(primaryIdx < 0 || (droppedBySupplementary != null && droppedBySupplementary[primaryIdx]))
        {
            return;
        }

        LiftBackResult primaryRes = resolved[primaryIdx];
        if(primaryRes == null || primaryRes.finalCigar() == null
                || primaryRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
            return;

        // The discriminator may have swapped the primary to an opposite-strand alt. The gate walks the read bases
        // against the genome in the placement's orientation, so reverse-complement them when the resolved strand
        // differs from the record's own (matches the per-candidate handling in reconcileAlignmentsToGenome).
        byte[] readBases = primary.getReadBases();
        if(readBases != null && primaryRes.negativeStrand() != primary.getReadNegativeStrandFlag())
        {
            readBases = Arrays.copyOf(readBases, readBases.length);
            SequenceUtil.reverseComplement(readBases);
        }

        // Peel the chosen primary (catches a fabricated terminal micro-junction merge created in a fresh M N M),
        // then reclaim a surviving terminal softclip - but only for a tx-match, and only here (post-resolve) so
        // supplementary resolve saw the original clip and a genomic over-clip is left untouched. A clean primary is unchanged.
        if(mOverhangGate != null && mOverhangGate.enabled())
        {
            OverhangGate.Result peeled = mOverhangGate.peel(
                    primaryRes.finalChrom(), primaryRes.finalPos(), primaryRes.finalCigar(), readBases);
            int newPos = peeled.pos();
            String newCigar = peeled.cigar();

            if(isTxMatchPrimary(primaryRes))
            {
                OverhangGate.Result reclaimed = mOverhangGate.reclaimTerminalSoftClip(
                        primaryRes.finalChrom(), newPos, newCigar, readBases);
                newPos = reclaimed.pos();
                newCigar = reclaimed.cigar();
            }

            if(newPos != primaryRes.finalPos() || !newCigar.equals(primaryRes.finalCigar()))
            {
                TARS_LOGGER.debug2("overhang-gate primary {}: {} -> {}",
                        primary.getReadName(), primaryRes.finalCigar(), newCigar);
                primaryRes = primaryRes.withRevisedCigar(
                        newPos, newCigar, newCigar.indexOf('N') >= 0,
                        primaryRes.updatedMapq(), "overhang-gated");
                resolved[primaryIdx] = primaryRes;
            }

            // Reclaim terminal softclips on tx-match XA alts as well. Supplementaries only ever merge into the
            // primary, so an alt's clip is never a splice the resolver still needs and is safe to consume here.
            // Genomic alts keep bwa's clip. Mutates the shared liftedAlignments list that buildLiftedXa emits from.
            reclaimTxMatchAlts(primaryRes.liftedAlignments(), primary);
        }
    }

    // Reclaim a surviving terminal softclip on each tx-match XA alt, mirroring the chosen primary's reclaim above.
    // Genomic alts (bwa's clip is a genuine mismatch call) and the primary are skipped. Read bases are taken in each
    // alt's placement orientation, matching the per-candidate handling in reconcileAlignmentsToGenome.
    private void reclaimTxMatchAlts(final List<LiftedAlignment> alignments, final SAMRecord primary)
    {
        byte[] forwardBases = primary.getReadBases();
        if(forwardBases == null || forwardBases.length == 0)
        {
            return;
        }

        boolean recordForward = !primary.getReadNegativeStrandFlag();
        byte[] reverseBases = null;

        for(int i = 0; i < alignments.size(); ++i)
        {
            LiftedAlignment alt = alignments.get(i);
            if(alt.IsPrimaryChoice || alt.Dropped || !alt.fromTxContig())
            {
                continue;
            }
            if(alt.LiftedCigar == null || alt.LiftedCigar.indexOf('S') < 0)
            {
                continue;
            }

            byte[] bases;
            if(alt.ForwardStrand == recordForward)
            {
                bases = forwardBases;
            }
            else
            {
                if(reverseBases == null)
                {
                    reverseBases = Arrays.copyOf(forwardBases, forwardBases.length);
                    SequenceUtil.reverseComplement(reverseBases);
                }
                bases = reverseBases;
            }

            OverhangGate.Result reclaimed = mOverhangGate.reclaimTerminalSoftClip(
                    alt.LiftedChrom, alt.LiftedPos, alt.LiftedCigar, bases);
            if(reclaimed.pos() != alt.LiftedPos || !reclaimed.cigar().equals(alt.LiftedCigar))
            {
                alignments.set(i, alt.withLiftedCigar(reclaimed.pos(), reclaimed.cigar()));
            }
        }
    }

    // Dedup key for supplementaries: lifted (chrom, pos, cigar) + lifted strand. Opposite-strand
    // placements at the same coords/cigar are kept as distinct records.
    private static String dedupKey(final LiftBackResult result, final SAMRecord record)
    {
        return result.finalChrom() + ":" + result.finalPos() + ":" + result.finalCigar()
                + ":" + (result.negativeStrand() ? '-' : '+');
    }

    // true if the result's lifted genomic span overlaps an excluded region. False for unmapped / unlifted results.
    private boolean liftsIntoExcludedRegion(final LiftBackResult result)
    {
        if(mExcludedRegions == null || result == null)
        {
            return false;
        }
        String cigar = result.finalCigar();
        if(result.finalChrom() == null || cigar == null || cigar.equals(SAMRecord.NO_ALIGNMENT_CIGAR))
        {
            return false;
        }
        int end = result.finalPos() + TextCigarCodec.decode(cigar).getReferenceLength() - 1;
        return mExcludedRegions.excludes(result.finalChrom(), result.finalPos(), end);
    }

    // SA entry key (chrom:pos:strand:cigar) of a supplementary's lifted placement, matching SaTagRewriter's
    // entry key so a dropped supp's entry can be removed from the primary's SA tag.
    private static String suppSaKey(final LiftBackResult result)
    {
        return result.finalChrom() + ":" + result.finalPos()
                + ":" + (result.negativeStrand() ? '-' : '+') + ":" + result.finalCigar();
    }

    private void applyAndWriteRecord(
            final SAMRecord record, final LiftBackResult result, final int nh, final boolean primaryPostProcessed,
            final Set<String> droppedSuppSaKeys, final LiftedMateInfoCache liftedMateInfoCache, final EmitSink sink)
    {
        mStats.record(record, result);

        LiftBackRecordOps.applyResultToRecord(record, result, liftedMateInfoCache);

        // A primary flipped to UNMAPPED post-lift (e.g. excluded region) still carries its mapped coords -
        // applyResultToRecord's UNMAPPED case is a no-op for input-unmapped reads, so unmap it explicitly.
        if(result.recordState() == RecordState.UNMAPPED && !record.getReadUnmappedFlag())
        {
            markPrimaryUnmapped(record);
        }

        String rewrittenSa = rewriteSaTag(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE), mResolver, droppedSuppSaKeys);

        // A record that stays supplementary must carry an SA tag -- REDUX dedup (FragmentCoords.fromRead)
        // reads the primary's coords from it. If every SA entry failed to lift the supplementary is orphaned
        // in genomic space and cannot be represented, so drop it rather than emit a supp with a null SA.
        if(rewrittenSa == null && record.getSupplementaryAlignmentFlag())
        {
            mStats.recordOrphanSuppDropped();
            return;
        }

        record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, rewrittenSa);
        patchMateFields(record, liftedMateInfoCache);

        if(result.recordState() != RecordState.UNMAPPED)
        {
            record.setAttribute("NH", nh);
        }

        // Over bwa's XA reporting cap: a primary that bwa emitted MAPQ 0 with no XA maps to too many loci to
        // place (the alt list was suppressed, not truncated). Unmap it. Keyed on the input state because the
        // missing XA makes the resolver see a single locus, which would otherwise bump MAPQ to 60.
        if(LiftBackRecordOps.exceedsMappingCap(result)
                && !record.getReadUnmappedFlag()
                && !record.getSupplementaryAlignmentFlag())
        {
            markPrimaryUnmapped(record);
            ++mOverCapUnmapped;
            TARS_LOGGER.debug2("over-cap unmap {}: inputMapq=0, no XA (past bwa XA cap)", record.getReadName());
        }

        // residual short-anchor primary: bwa scored it below the default -T 30 floor and liftback couldn't
        // improve it (not resolved / extended / collapsed). Unmap it so the noise alignment doesn't pollute
        // the locus. Gated on supplementary resolve running, matching the supp drop's -T 19 rationale.
        if(mSupplementaryResolver != null && !primaryPostProcessed
                && !record.getReadUnmappedFlag()
                && !record.getSupplementaryAlignmentFlag())
        {
            Integer alignmentScore = record.getIntegerAttribute(AS_TAG);
            if(alignmentScore != null && alignmentScore < PRIMARY_AS_UNMAP_THRESHOLD)
            {
                markPrimaryUnmapped(record);
                mStats.recordLowAsPrimaryUnmapped();
                TARS_LOGGER.debug2("AS-floor unmap {}: AS={} < {}", record.getReadName(), alignmentScore, PRIMARY_AS_UNMAP_THRESHOLD);
            }
        }

        refreshNmDropMd(record, mRefSource);

        // A placement TARS made confident on its own (bumped a bwa MAPQ-0 read, or swapped to a different locus) but
        // whose recomputed genomic NM is a large fraction of the read length is not a trustworthy unique placement
        // (e.g. a large-D lift artifact heavily mismatching the genome) -- revert it to MAPQ 0 rather than assert
        // confidence. bwa's own kept confident calls (input MAPQ > 0, not swapped) are never touched.
        if((result.inputMapq() == 0 || result.swapped())
                && !record.getReadUnmappedFlag() && !record.getSupplementaryAlignmentFlag()
                && record.getMappingQuality() >= CONFIDENT_MAPQ)
        {
            Integer editDistance = record.getIntegerAttribute("NM");
            byte[] readBases = record.getReadBases();
            if(editDistance != null && readBases != null && readBases.length > 0
                    && editDistance > MAX_CONFIDENT_NM_FRACTION * readBases.length)
            {
                record.setMappingQuality(0);
                mStats.recordLowIdentityDemoted();
                TARS_LOGGER.debug2("low-identity demote {}: NM={} over {}% of {}bp",
                        record.getReadName(), editDistance, (int) (MAX_CONFIDENT_NM_FRACTION * 100), readBases.length);
            }
        }

        sink.emit(record, result);
    }
}
