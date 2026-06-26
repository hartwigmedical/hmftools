package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.LiftBackRecordOps.markPrimaryUnmapped;
import static com.hartwig.hmftools.tars.liftback.LiftBackRecordOps.refreshNmDropMd;
import static com.hartwig.hmftools.tars.liftback.MateFieldPatcher.patchMateFields;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.rewriteSaTag;
import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.SUPP_AS_DROP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.SUPP_RESCUE_MAPQ_CAP;
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
import com.hartwig.hmftools.tars.liftback.rescue.JunctionRescueResolver;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.rescue.RescueCandidate;
import com.hartwig.hmftools.tars.liftback.rescue.RescueResult;
import com.hartwig.hmftools.tars.liftback.rescue.RescueStatistics;
import com.hartwig.hmftools.tars.liftback.rescue.RescueSupplementary;
import com.hartwig.hmftools.tars.liftback.tailextend.TailExtensionStatistics;
import com.hartwig.hmftools.tars.liftback.tailextend.TerminalReconciler;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;

// Stateless-per-group transform engine: resolves one contiguous read-name group to genomic coordinates,
// runs the optional rescue / terminal-collapse / tail-extend passes, patches mate fields against a
// provided LiftedMateInfoCache, and emits each kept record to an EmitSink. The cache is supplied by the
// caller so the standalone TarsApplication can pass its whole-sample pass-1 cache while the REDUX worker
// passes a fresh per-group cache (a complete name-sorted group is self-sufficient -- see
// [[project_redux_migration_scoping]]). Extracted verbatim from TarsApplication so both paths share one
// implementation.
public class LiftBackGroupProcessor
{
    private final LiftBackResolver mResolver;

    // optional rescue resolver: merges primary + supp across annotated junctions. Null when disabled.
    private final JunctionRescueResolver mJunctionRescueResolver;

    // null when extend-softclip-tails / terminal-collapse are disabled. Holds both terminal passes.
    private final TerminalReconciler mTerminalReconciler;
    private final JunctionCanonicalizer mJunctionCanonicalizer;

    // genomic reference for the post-lift NM recompute. Null when no ref is loaded (rescue + extend both
    // off), in which case NM is cleared rather than recomputed.
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
            final LiftBackResolver resolver, final JunctionRescueResolver junctionRescueResolver,
            final TerminalReconciler terminalReconciler,
            final JunctionCanonicalizer junctionCanonicalizer, final RefSequenceSource refSource,
            final ExcludedRegions excludedRegions, final LiftBackStats stats)
    {
        mResolver = resolver;
        mJunctionRescueResolver = junctionRescueResolver;
        mTerminalReconciler = terminalReconciler;
        mJunctionCanonicalizer = junctionCanonicalizer;
        mRefSource = refSource;
        mExcludedRegions = excludedRegions;
        mStats = stats;
    }

    // Processes all records sharing one read name as a group so primary + supps resolve together. This
    // lets rescue see all split-read components and lets /2 hint /1's introns (and vice versa).
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

        // Two passes per pair so mate /2 can use mate /1's rescued junctions as a hint, and a
        // re-decide of mate /1 picks up mate /2's hints when mate /1 wasn't initially rescued.
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
            // pass-effect counters it bumped (collapse / tail-extend / canonicalize / rescue / excluded reads) so
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

    // Push the post-rescue primary back into the mate cache so the partner mate's MC tag uses
    // the merged/extended cigar instead of the pre-rescue one the seed used.
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
        cache.recordPrimaryAlignment(primary.getReadName(), primary.getFirstOfPairFlag(), refreshed);
    }

    // Snapshot of the pass-effect counters mutated by one decideMateGroup pass, so a discarded provisional mate
    // decision can be rolled back rather than double-counted in the summary (see processNameGroup).
    private static final class PassCounterSnapshot
    {
        final int[] Tail;
        final int[] Rescue;
        final long CollapseLeading;
        final long CollapseTrailing;
        final long Canonicalized;
        final long ExcludedReads;

        PassCounterSnapshot(
                final int[] tail, final int[] rescue, final long collapseLeading, final long collapseTrailing,
                final long canonicalized, final long excludedReads)
        {
            Tail = tail;
            Rescue = rescue;
            CollapseLeading = collapseLeading;
            CollapseTrailing = collapseTrailing;
            Canonicalized = canonicalized;
            ExcludedReads = excludedReads;
        }
    }

    private PassCounterSnapshot snapshotPassCounters()
    {
        return new PassCounterSnapshot(
                mTerminalReconciler != null ? mTerminalReconciler.statistics().snapshot() : null,
                mJunctionRescueResolver != null ? mJunctionRescueResolver.statistics().snapshot() : null,
                mTerminalReconciler != null ? mTerminalReconciler.collapsedLeading() : 0,
                mTerminalReconciler != null ? mTerminalReconciler.collapsedTrailing() : 0,
                mJunctionCanonicalizer != null ? mJunctionCanonicalizer.junctionsShifted() : 0,
                mExcludedReads);
    }

    // Subtract the (after - before) delta a discarded provisional pass contributed, leaving only the kept counts.
    private void rewindProvisionalCounters(final PassCounterSnapshot before, final PassCounterSnapshot after)
    {
        if(mTerminalReconciler != null)
        {
            TailExtensionStatistics tail = mTerminalReconciler.statistics();
            int[] current = tail.snapshot();
            for(int i = 0; i < current.length; ++i)
            {
                current[i] -= after.Tail[i] - before.Tail[i];
            }
            tail.restore(current);
            mTerminalReconciler.restoreCollapseCounters(
                    mTerminalReconciler.collapsedLeading() - (after.CollapseLeading - before.CollapseLeading),
                    mTerminalReconciler.collapsedTrailing() - (after.CollapseTrailing - before.CollapseTrailing));
        }
        if(mJunctionRescueResolver != null)
        {
            RescueStatistics rescue = mJunctionRescueResolver.statistics();
            int[] current = rescue.snapshot();
            for(int i = 0; i < current.length; ++i)
            {
                current[i] -= after.Rescue[i] - before.Rescue[i];
            }
            rescue.restore(current);
        }
        if(mJunctionCanonicalizer != null)
        {
            mJunctionCanonicalizer.restoreJunctionsShifted(
                    mJunctionCanonicalizer.junctionsShifted() - (after.Canonicalized - before.Canonicalized));
        }
        mExcludedReads -= after.ExcludedReads - before.ExcludedReads;
    }

    // Per-mate decision: resolves every record, runs the rescue resolver, returns the lifted results +
    // per-record drop flags + introns introduced by rescue (for mate hinting). Does not emit.
    private static final class MateDecision
    {
        LiftBackResult[] Resolved;
        boolean[] DroppedByRescue;
        LiftBackResult PrimaryResult;
        List<ChrBaseRegion> IntroducedIntrons;
        int PrimaryIndex;
        boolean PrimaryPostProcessed;

        MateDecision(final LiftBackResult[] resolved, final boolean[] droppedByRescue,
                final LiftBackResult primaryResult, final List<ChrBaseRegion> introducedIntrons,
                final int primaryIndex, final boolean primaryPostProcessed)
        {
            Resolved = resolved;
            DroppedByRescue = droppedByRescue;
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

        // a valid BAM always has a primary in the mate group, so resolve it directly. Candidate cigars are
        // reconciled (collapse + tail-extend) before the discriminator runs - see reconcileAlignmentsToGenome.
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
        boolean[] droppedByRescue = applyJunctionRescue(records, resolved, primary, mateHintIntrons, introduced);

        // Finalize the chosen primary: reconcile its lifted cigar (terminal collapse + tail-extend), then
        // canonicalize. Not a redundant repeat of the per-candidate reconciliation in reconcileAlignmentsToGenome -
        // it does real work in two cases:
        //  1. rescue merged a primary+supp into a fresh M N M cigar that never existed at candidate time, so it
        //     can carry its own fabricated terminal micro-junction or reclaimable softclip; and
        //  2. the discriminator swapped to a cross-chromosome alt, which reconcileAlignmentsToGenome collapsed but
        //     did NOT tail-extend (tail-extend is scoped to self + co-located alts) - this is its first tail-extend.
        // Case 2 is why that scoping is output-safe: whatever primary wins is fully reconciled here regardless.
        // A primary already fully reconciled at candidate time skips every pass via the N-present / S-present guards.
        reconcileChosenPrimary(records, resolved, primary, droppedByRescue);

        // post-lift exclusion: if the final primary placement lands in an excluded (rRNA / contamination) region,
        // unmap it REDUX-style - flip the result to UNMAPPED so the mate is coordinated via the cache (willBeUnmapped)
        // and the record is unmapped at apply time. Lifted genomic coords are the only space tx-contig reads can be
        // tested against the genomic region list, so this is post-lift, not the old pre-lift fragment pre-filter.
        if(primaryIdx >= 0 && liftsIntoExcludedRegion(resolved[primaryIdx]))
        {
            resolved[primaryIdx] = LiftBackResult.unmapped(LiftBackResult.RecordRole.PRIMARY, "excluded_region_unmapped");
            ++mExcludedReads;
        }

        // rescue / collapse / tail-extend each replace resolved[primaryIdx] with a new result object when
        // they improve the primary. A changed reference means liftback post-processed it, so its stale bwa
        // AS must not be used to AS-filter it (see PRIMARY_AS_UNMAP_THRESHOLD).
        boolean primaryPostProcessed = primaryIdx >= 0 && resolved[primaryIdx] != primaryResult;

        return new MateDecision(resolved, droppedByRescue,
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
        boolean[] droppedByRescue = decision.DroppedByRescue;
        SAMRecord primary = null;
        for(SAMRecord r : records)
            if(!r.getSupplementaryAlignmentFlag())
            {
                primary = r;
                break;
            }

        // Dedup supplementaries that lift to the same (chrom, pos, cigar, strand) -- bwa can emit the
        // same junction across multiple transcript contigs and they collapse after liftback.
        Set<String> emittedSuppKeys = new HashSet<>();
        boolean[] willEmit = new boolean[records.size()];
        for(int i = 0; i < records.size(); i++)
        {
            SAMRecord record = records.get(i);
            LiftBackResult result = resolved[i];
            boolean drop = droppedByRescue != null && droppedByRescue[i];
            if(!drop && record != primary && record.getSupplementaryAlignmentFlag())
            {
                String key = dedupKey(result, record);
                if(!emittedSuppKeys.add(key))
                {
                    drop = true;
                }
            }
            // Drop supps the rescue pass left behind that exist only because bwa-mem2 was run with -T 19
            // below its default of 30 (see SUPP_AS_DROP_THRESHOLD). Only applies when rescue ran, because
            // a configuration without rescue might want to retain these supps for other reasons.
            if(!drop && mJunctionRescueResolver != null && record.getSupplementaryAlignmentFlag()
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

    // builds a RescueCandidate from the post-lift primary + its supplementary records and applies the
    // merge if accepted. Returns a parallel array marking which records were absorbed by the merge and
    // should be dropped from emission. Returns null when rescue is disabled / no-op. introducedIntronsOut
    // is appended-to (caller-provided list) so the partner mate can use them as hints in the second pass.
    private boolean[] applyJunctionRescue(
            final List<SAMRecord> records, final LiftBackResult[] resolved, final SAMRecord primary,
            final List<ChrBaseRegion> mateHintIntrons, final List<ChrBaseRegion> introducedIntronsOut)
    {
        if(mJunctionRescueResolver == null || primary == null || primary.getReadUnmappedFlag())
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

        List<RescueSupplementary> suppDtos = new ArrayList<>(suppIndices.size());
        for(int i = 0; i < suppIndices.size(); i++)
        {
            int idx = suppIndices.get(i);
            SAMRecord supp = records.get(idx);
            LiftBackResult suppRes = resolved[idx];
            if(suppRes == null || suppRes.finalCigar() == null
                    || suppRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
                continue;
            suppDtos.add(new RescueSupplementary(
                    i, suppRes.finalChrom(), !supp.getReadNegativeStrandFlag(),
                    suppRes.finalPos(), suppRes.finalCigar(), supp.getMappingQuality()));
        }

        RescueCandidate candidate = new RescueCandidate(
                primaryRes.finalChrom(), !primary.getReadNegativeStrandFlag(), primary.getReadLength(),
                primaryRes.finalPos(), primaryRes.finalCigar(), primary.getMappingQuality(), suppDtos,
                primary.getReadBases(), mateHintIntrons);

        RescueResult result = mJunctionRescueResolver.resolve(candidate);
        if(!result.merged())
        {
            return null;
        }

        // rewrite the primary's LiftBackResult with the merged (now-spliced) cigar + start; the start may
        // shift for a left-extend, and MAPQ caps at SUPP_RESCUE_MAPQ_CAP to mark it a constructed alignment.
        // mark merged supps for drop.
        resolved[primaryIdx] = primaryRes.withRevisedCigar(
                result.mergedStart(), result.mergedCigar(), true,
                Math.min(primaryRes.updatedMapq(), SUPP_RESCUE_MAPQ_CAP), "rescued-via-supp");
        TARS_LOGGER.debug2("rescue-merged {}: {}:{} {} -> {} (depth {})",
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

    // Reconciles candidate cigars against the genome before the discriminator: collapse a fabricated terminal
    // micro-junction (sub-trust anchor across an N), then walk a reclaimable terminal softclip into contiguous
    // genome. This turns the discriminator's RefSoftClipped / RefFullMatch / TxHasNCigar inputs from raw-bwa
    // guesses into measured facts. Read bases are taken in each candidate's own orientation; self and same-strand
    // alts use the record bases as-is, opposite-strand alts use the reverse complement.
    //
    // Collapse runs on every candidate: it can shift a start by a whole intron, merging a fabricated far
    // placement back onto self's locus (turns a spurious multimapper into a single-locus rescue). Tail-extend
    // is restricted to co-located alts (NOT self): it shifts a start by at most a read's tail, so it never merges
    // loci, and it only feeds the single-locus RefSoftClipped/RefFullMatch contest - running it on far-locus
    // cross-locus alts is wasted ref-genome work (the dominant cost on multimapper-heavy samples). Self is left
    // collapse-only so its terminal softclip survives for junction rescue (a primary+supp merge needs the clip)
    // and for reconcileChosenPrimary's post-rescue tail-extend - i.e. self's tail-extend runs after rescue, as
    // it did before this pass was hoisted ahead of the discriminator.
    private void reconcileAlignmentsToGenome(final List<LiftedAlignment> alignments, final SAMRecord record)
    {
        if(mTerminalReconciler == null || alignments.isEmpty())
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

        // Fast path for the dominant single-candidate read: reconcile self directly and skip the per-placement
        // dedup map. Self is collapse-only (allowTailExtend false) - its softclip is left for rescue and the
        // post-rescue tail-extend in reconcileChosenPrimary. Self uses the record's own orientation (forward bases).
        if(alignments.size() == 1)
        {
            if(self.LiftedCigar != null)
            {
                ReconciledPlacement placement = reconcileCigarToGenome(self, forwardBases, false);
                if(placement.Pos != self.LiftedPos || !placement.Cigar.equals(self.LiftedCigar))
                {
                    alignments.set(0, self.withLiftedCigar(placement.Pos, placement.Cigar));
                }
            }
            return;
        }

        byte[] reverseBases = null;

        // Many candidates lift to an identical genomic placement (e.g. the packed isoform contigs of one gene all
        // collapsing to one locus). The collapse / tail-extend ref-genome lookups are deterministic in
        // (chrom, pos, cigar, strand) -- allowTailExtend is just chrom == self's chrom, so it is constant per key --
        // so each distinct placement is reconciled once and the result reused for every duplicate.
        Map<String, ReconciledPlacement> placementsByKey = new HashMap<>();

        for(int i = 0; i < alignments.size(); ++i)
        {
            LiftedAlignment alignment = alignments.get(i);
            if(alignment.LiftedCigar == null)
            {
                continue;
            }

            String key = alignment.LiftedChrom + ':' + alignment.LiftedPos + ':' + alignment.LiftedCigar
                    + (alignment.ForwardStrand ? '+' : '-');
            ReconciledPlacement placement = placementsByKey.get(key);
            if(placement == null)
            {
                boolean allowTailExtend = alignment != self && coLocated(alignment, self);

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

                placement = reconcileCigarToGenome(alignment, bases, allowTailExtend);
                placementsByKey.put(key, placement);
            }

            if(placement.Pos != alignment.LiftedPos || !placement.Cigar.equals(alignment.LiftedCigar))
            {
                alignments.set(i, alignment.withLiftedCigar(placement.Pos, placement.Cigar));
            }
        }
    }

    // Collapse a fabricated terminal micro-junction, then walk a reclaimable terminal softclip into contiguous
    // genome. Reads only (chrom, pos, cigar, bases, allowTailExtend), so the result is cacheable per placement.
    private ReconciledPlacement reconcileCigarToGenome(
            final LiftedAlignment alignment, final byte[] bases, final boolean allowTailExtend)
    {
        TerminalReconciler.ReconcileResult result =
                mTerminalReconciler.reconcile(alignment.LiftedChrom, alignment.LiftedPos, alignment.LiftedCigar, bases, allowTailExtend);
        return new ReconciledPlacement(result.pos(), result.cigar());
    }

    // Cached collapse/tail-extend output for one lifted placement (see reconcileAlignmentsToGenome).
    private static final class ReconciledPlacement
    {
        final int Pos;
        final String Cigar;

        ReconciledPlacement(final int pos, final String cigar)
        {
            Pos = pos;
            Cigar = cigar;
        }
    }

    // Same chromosome as self: tail-extend reaches any alt that could share self's locus after a short walk,
    // while still skipping the cross-chromosome alts of a multimapper (the bulk of the wasted ref
    // lookups). A pos window proved too tight - it dropped same-chrom contest alts and lost their rescues.
    private static boolean coLocated(final LiftedAlignment alignment, final LiftedAlignment self)
    {
        return alignment.LiftedChrom.equals(self.LiftedChrom);
    }

    // Finalize the chosen primary in one pass: terminal micro-junction collapse, then softclip tail-extend, then
    // junction canonicalize, each operating on the running result. Runs on whatever primary the discriminator or
    // rescue chose (see the call site for why it is not rescue-specific). A clean primary skips each pass via the
    // N-present / S-present guards, so it costs only the lookup. Same sequence/gating as the three former passes.
    private void reconcileChosenPrimary(
            final List<SAMRecord> records, final LiftBackResult[] resolved, final SAMRecord primary,
            final boolean[] droppedByRescue)
    {
        if(mTerminalReconciler == null && mJunctionCanonicalizer == null)
        {
            return;
        }
        if(primary == null || primary.getReadUnmappedFlag())
        {
            return;
        }

        int primaryIdx = indexOfPrimary(records, primary);
        if(primaryIdx < 0 || (droppedByRescue != null && droppedByRescue[primaryIdx]))
        {
            return;
        }

        LiftBackResult primaryRes = resolved[primaryIdx];
        if(primaryRes == null || primaryRes.finalCigar() == null
                || primaryRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
            return;

        // The discriminator may have swapped the primary to an opposite-strand alt. tryCollapse / tryExtend walk
        // the read bases against the genome in the placement's orientation, so reverse-complement them when the
        // resolved strand differs from the record's own (matches the per-candidate handling in reconcileAlignmentsToGenome).
        byte[] readBases = primary.getReadBases();
        if(readBases != null && primaryRes.negativeStrand() != primary.getReadNegativeStrandFlag())
        {
            readBases = Arrays.copyOf(readBases, readBases.length);
            SequenceUtil.reverseComplement(readBases);
        }

        // 1 + 2. reconcile the chosen primary: terminal micro-junction collapse, then softclip tail-extend, via
        // the single TerminalReconciler.reconcile() entry point (one definition shared with the per-candidate
        // pass, so the two cannot drift). This is the chosen primary's tail-extend, deferred from candidate time
        // so rescue saw the original softclip (see reconcileAlignmentsToGenome). A clean primary is unchanged.
        if(mTerminalReconciler != null)
        {
            TerminalReconciler.ReconcileResult reconciled = mTerminalReconciler.reconcile(
                    primaryRes.finalChrom(), primaryRes.finalPos(), primaryRes.finalCigar(), readBases, true);
            if(reconciled.pos() != primaryRes.finalPos() || !reconciled.cigar().equals(primaryRes.finalCigar()))
            {
                TARS_LOGGER.debug2("reconcile-primary {}: {} -> {}",
                        primary.getReadName(), primaryRes.finalCigar(), reconciled.cigar());
                primaryRes = primaryRes.withRevisedCigar(
                        reconciled.pos(), reconciled.cigar(), reconciled.cigar().indexOf('N') >= 0,
                        primaryRes.updatedMapq(), "reconciled");
                resolved[primaryIdx] = primaryRes;
            }
        }

        // 3. slide an interior junction onto a canonical splice motif (N-cigars only)
        if(mJunctionCanonicalizer != null && primaryRes.hasNCigar())
        {
            JunctionCanonicalizationResult canon = mJunctionCanonicalizer.tryCanonicalize(
                    primaryRes.finalChrom(), primaryRes.finalPos(), primaryRes.finalCigar(), readBases);
            if(canon.changed())
            {
                TARS_LOGGER.debug2("canonicalize {}: {} -> {}",
                        primary.getReadName(), primaryRes.finalCigar(), canon.newCigar());
                primaryRes = primaryRes.withRevisedCigar(
                        primaryRes.finalPos(), canon.newCigar(), primaryRes.hasNCigar(),
                        primaryRes.updatedMapq(), "junction-canonicalized");
                resolved[primaryIdx] = primaryRes;
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
        // missing XA makes the resolver see a single locus, which would otherwise rescue MAPQ to 60.
        if(LiftBackRecordOps.exceedsMappingCap(result)
                && !record.getReadUnmappedFlag()
                && !record.getSupplementaryAlignmentFlag())
        {
            markPrimaryUnmapped(record);
            ++mOverCapUnmapped;
            TARS_LOGGER.debug2("over-cap unmap {}: inputMapq=0, no XA (past bwa XA cap)", record.getReadName());
        }

        // residual short-anchor primary: bwa scored it below the default -T 30 floor and liftback couldn't
        // improve it (not rescued / extended / collapsed). Unmap it so the noise alignment doesn't pollute
        // the locus. Gated on rescue running, matching the supp drop's -T 19 rationale.
        if(mJunctionRescueResolver != null && !primaryPostProcessed
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
        sink.emit(record, result);
    }
}
