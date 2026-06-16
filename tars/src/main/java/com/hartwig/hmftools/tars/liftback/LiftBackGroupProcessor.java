package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.LiftBackRecordOps.markPrimaryUnmapped;
import static com.hartwig.hmftools.tars.liftback.LiftBackRecordOps.refreshNmDropMd;
import static com.hartwig.hmftools.tars.liftback.MateFieldPatcher.patchMateFields;
import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.SA_ATTRIBUTE;
import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.rewriteSaTag;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.rescue.JunctionRescueResolver;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.rescue.RescueCandidate;
import com.hartwig.hmftools.tars.liftback.rescue.RescueResult;
import com.hartwig.hmftools.tars.liftback.rescue.RescueSupplementary;
import com.hartwig.hmftools.tars.liftback.tailextend.SoftclipTailExtender;
import com.hartwig.hmftools.tars.liftback.tailextend.TailExtensionResult;
import com.hartwig.hmftools.tars.liftback.tailextend.TerminalCollapseResult;
import com.hartwig.hmftools.tars.liftback.tailextend.TerminalMicroJunctionCollapser;

import htsjdk.samtools.SAMRecord;

// Stateless-per-group transform engine: resolves one contiguous read-name group to genomic coordinates,
// runs the optional rescue / terminal-collapse / tail-extend passes, patches mate fields against a
// provided LiftedMateInfoCache, and emits each kept record to an EmitSink. The cache is supplied by the
// caller so the standalone SpliceLiftBack can pass its whole-sample pass-1 cache while the REDUX worker
// passes a fresh per-group cache (a complete name-sorted group is self-sufficient -- see
// [[project_redux_migration_scoping]]). Extracted verbatim from SpliceLiftBack so both paths share one
// implementation.
public class LiftBackGroupProcessor
{
    private final LiftBackResolver mResolver;

    // optional rescue resolver: merges primary + supp across annotated junctions. Null when disabled.
    private final JunctionRescueResolver mJunctionRescueResolver;

    // null when extend-softclip-tails / terminal-collapse are disabled.
    private final SoftclipTailExtender mSoftclipTailExtender;
    private final TerminalMicroJunctionCollapser mTerminalJunctionCollapser;
    private final JunctionCanonicalizer mJunctionCanonicalizer;

    // genomic reference for the post-lift NM recompute. Null when no ref is loaded (rescue + extend both
    // off), in which case NM is cleared rather than recomputed.
    private final RefSequenceSource mRefSource;

    private final int mUnmapAboveNh;
    private final int mUnmapBelowMapq;

    private final LiftBackStats mStats;

    private static final String AS_TAG = "AS";

    // bwa-mem2's default -T (minimum alignment score) is 30. The RNA splice run uses -T 19 deliberately
    // to surface short-anchor supplementary records that JunctionRescueResolver can merge across annotated
    // junctions. Supps in the [19, 30) AS band that survive the rescue pass are residual noise -- drop
    // them from the output BAM after rescue + tail-extend have had their chance.
    static final int SUPP_AS_DROP_THRESHOLD = 30;

    // Same bwa -T 19 rationale for primaries: a primary whose AS is still below the default floor of 30
    // AFTER liftback (and which rescue/tail-extend/collapse did not improve) is a residual short-anchor
    // alignment. AS is never recomputed in liftback, so post-processed primaries keep a stale-low AS and
    // must be exempt. Such primaries are unmapped (not dropped) to keep the pair + SA references intact,
    // matching the unmap_below_mapq policy.
    static final int PRIMARY_AS_UNMAP_THRESHOLD = 30;

    // BWA emits MAPQ=60 for a clean unique alignment. Rescued primaries are constructed by us, not
    // directly emitted by BWA, so we cap at 55 to signal a primary+supp merge. Capping never goes UP.
    private static final int RESCUED_MAPQ_CAP = 55;

    // sink for emitted records: the standalone writes a BAM + TSV, the REDUX worker writes the shared BAM.
    public interface EmitSink
    {
        void emit(SAMRecord record, LiftBackResult result);
    }

    public LiftBackGroupProcessor(
            final LiftBackResolver resolver, final JunctionRescueResolver junctionRescueResolver,
            final SoftclipTailExtender softclipTailExtender, final TerminalMicroJunctionCollapser terminalJunctionCollapser,
            final JunctionCanonicalizer junctionCanonicalizer, final RefSequenceSource refSource,
            final int unmapAboveNh, final int unmapBelowMapq, final LiftBackStats stats)
    {
        mResolver = resolver;
        mJunctionRescueResolver = junctionRescueResolver;
        mSoftclipTailExtender = softclipTailExtender;
        mTerminalJunctionCollapser = terminalJunctionCollapser;
        mJunctionCanonicalizer = junctionCanonicalizer;
        mRefSource = refSource;
        mUnmapAboveNh = unmapAboveNh;
        mUnmapBelowMapq = unmapBelowMapq;
        mStats = stats;
    }

    // Processes all records sharing one read name as a group so primary + supps resolve together. This
    // lets rescue see all split-read components and lets /2 hint /1's introns (and vice versa).
    public void processNameGroup(
            final List<SAMRecord> group, final LiftedMateInfoCache liftedMateInfoCache, final EmitSink sink)
    {
        final List<SAMRecord> firstOfPair = new ArrayList<>();
        final List<SAMRecord> secondOfPair = new ArrayList<>();
        for(final SAMRecord record : group)
        {
            if(!record.getReadPairedFlag() || record.getFirstOfPairFlag())
                firstOfPair.add(record);
            else
                secondOfPair.add(record);
        }

        // Two passes per pair so mate /2 can use mate /1's rescued junctions as a hint, and a
        // re-decide of mate /1 picks up mate /2's hints when mate /1 wasn't initially rescued.
        final MateDecision m1 = decideMateGroup(firstOfPair, Collections.emptyList());
        refreshMateInfoCache(firstOfPair, m1, liftedMateInfoCache);
        final MateDecision m2 = decideMateGroup(secondOfPair, m1.IntroducedIntrons);
        refreshMateInfoCache(secondOfPair, m2, liftedMateInfoCache);
        final MateDecision m1Final = (m1.IntroducedIntrons.isEmpty() && !m2.IntroducedIntrons.isEmpty())
                ? decideMateGroup(firstOfPair, m2.IntroducedIntrons)
                : m1;
        if(m1Final != m1)
            refreshMateInfoCache(firstOfPair, m1Final, liftedMateInfoCache);

        writeMateGroup(firstOfPair, m1Final, liftedMateInfoCache, sink);
        writeMateGroup(secondOfPair, m2, liftedMateInfoCache, sink);
    }

    // Push the post-rescue primary back into the mate cache so the partner mate's MC tag uses
    // the merged/extended cigar instead of the pre-rescue one the seed used.
    private void refreshMateInfoCache(
            final List<SAMRecord> mateRecords, final MateDecision decision, final LiftedMateInfoCache cache)
    {
        if(mateRecords.isEmpty() || decision.PrimaryResult == null)
            return;
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
            return;
        final LiftedMateInfo refreshed = LiftBackRecordOps.toLiftedMateInfo(
                primary, decision.PrimaryResult, mUnmapAboveNh, mUnmapBelowMapq);
        cache.recordPrimaryAlignment(primary.getReadName(), primary.getFirstOfPairFlag(), refreshed);
    }

    // Per-mate decision: resolves every record, runs the rescue resolver, returns the lifted results +
    // per-record drop flags + introns introduced by rescue (for mate hinting). Does not emit.
    private static final class MateDecision
    {
        final LiftBackResult[] Resolved;
        final boolean[] DroppedByRescue;
        final LiftBackResult PrimaryResult;
        final List<ChrBaseRegion> IntroducedIntrons;
        final int PrimaryIndex;
        final boolean PrimaryPostProcessed;

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
            return MateDecision.empty();

        SAMRecord primary = null;
        for(final SAMRecord record : records)
        {
            if(record.getSupplementaryAlignmentFlag())
                continue;
            if(primary == null)
                primary = record;
        }

        final LiftBackResult primaryResult = primary != null ? mResolver.resolve(primary) : null;

        final LiftBackResult[] resolved = new LiftBackResult[records.size()];
        for(int i = 0; i < records.size(); i++)
        {
            final SAMRecord record = records.get(i);
            resolved[i] = (record == primary && primaryResult != null) ? primaryResult : mResolver.resolve(record);
        }

        final int primaryIdx = primary != null ? indexOfPrimary(records, primary) : -1;

        final List<ChrBaseRegion> introduced = new ArrayList<>();
        final boolean[] droppedByRescue = applyJunctionRescue(records, resolved, primary, mateHintIntrons, introduced);

        // Collapse a spurious tx-contig terminal micro-junction before tail extension, so the extender
        // sees the corrected (contiguous) cigar rather than the fabricated junction.
        applyTerminalJunctionCollapse(records, resolved, primary, droppedByRescue);

        applyTailExtension(records, resolved, primary, droppedByRescue);

        // Last lifted-cigar pass: slide any interior junction left a few bp onto a canonical splice motif
        // (fixes a tx-contig boundary deletion that lifted the intron off the true GT-AG site).
        applyJunctionCanonicalization(records, resolved, primary, droppedByRescue);

        // rescue / collapse / tail-extend each replace resolved[primaryIdx] with a new result object when
        // they improve the primary. A changed reference means liftback post-processed it, so its stale bwa
        // AS must not be used to AS-filter it (see PRIMARY_AS_UNMAP_THRESHOLD).
        final boolean primaryPostProcessed = primaryIdx >= 0 && resolved[primaryIdx] != primaryResult;

        return new MateDecision(resolved, droppedByRescue,
                primaryIdx >= 0 ? resolved[primaryIdx] : primaryResult,
                introduced, primaryIdx, primaryPostProcessed);
    }

    private static int indexOfPrimary(final List<SAMRecord> records, final SAMRecord primary)
    {
        for(int i = 0; i < records.size(); ++i)
            if(records.get(i) == primary)
                return i;
        return -1;
    }

    private void writeMateGroup(
            final List<SAMRecord> records, final MateDecision decision,
            final LiftedMateInfoCache liftedMateInfoCache, final EmitSink sink)
    {
        if(records.isEmpty())
            return;

        final LiftBackResult[] resolved = decision.Resolved;
        final boolean[] droppedByRescue = decision.DroppedByRescue;
        SAMRecord primary = null;
        for(SAMRecord r : records)
            if(!r.getSupplementaryAlignmentFlag()) { primary = r; break; }

        // Dedup supplementaries that lift to the same (chrom, pos, cigar, strand) -- bwa can emit the
        // same junction across multiple transcript contigs and they collapse after liftback.
        final Set<String> emittedSuppKeys = new HashSet<>();
        final boolean[] willEmit = new boolean[records.size()];
        for(int i = 0; i < records.size(); i++)
        {
            final SAMRecord record = records.get(i);
            final LiftBackResult result = resolved[i];
            boolean drop = droppedByRescue != null && droppedByRescue[i];
            if(!drop && record != primary && record.getSupplementaryAlignmentFlag())
            {
                final String key = dedupKey(result, record);
                if(!emittedSuppKeys.add(key))
                    drop = true;
            }
            // Drop supps the rescue pass left behind that exist only because bwa-mem2 was run with -T 19
            // below its default of 30 (see SUPP_AS_DROP_THRESHOLD). Only applies when rescue ran, because
            // a configuration without rescue might want to retain these supps for other reasons.
            if(!drop && mJunctionRescueResolver != null && record.getSupplementaryAlignmentFlag()
                    && !record.getReadUnmappedFlag())
            {
                final Integer alignmentScore = record.getIntegerAttribute(AS_TAG);
                if(alignmentScore != null && alignmentScore < SUPP_AS_DROP_THRESHOLD)
                {
                    drop = true;
                    mStats.recordLowAsSuppDropped();
                }
            }
            willEmit[i] = !drop;
        }

        // NH = the number of distinct genomic loci the read lifts back to, taken from the resolver's
        // chrom:pos-keyed locus count. Counting emitted records would inflate NH by tx-contig
        // representation multiplicity (one junction repeated across many transcript contigs all lifting
        // to the same locus).
        final int nh = decision.PrimaryResult != null ? Math.max(decision.PrimaryResult.numLoci(), 1) : 1;

        for(int i = 0; i < records.size(); i++)
        {
            final SAMRecord record = records.get(i);
            final LiftBackResult result = resolved[i];
            if(!willEmit[i])
            {
                mStats.record(record, result);
                continue;
            }
            final boolean primaryPostProcessed = i == decision.PrimaryIndex && decision.PrimaryPostProcessed;
            applyAndWriteRecord(record, result, nh, primaryPostProcessed, liftedMateInfoCache, sink);
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
            return null;

        int primaryIdx = -1;
        final List<Integer> suppIndices = new ArrayList<>();
        for(int i = 0; i < records.size(); i++)
        {
            final SAMRecord r = records.get(i);
            if(r == primary)
                primaryIdx = i;
            else if(r.getSupplementaryAlignmentFlag() && !r.getReadUnmappedFlag())
                suppIndices.add(i);
        }
        if(primaryIdx < 0)
            return null;

        final LiftBackResult primaryRes = resolved[primaryIdx];
        if(primaryRes == null || primaryRes.finalCigar() == null
                || primaryRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
            return null;

        final List<RescueSupplementary> suppDtos = new ArrayList<>(suppIndices.size());
        for(int i = 0; i < suppIndices.size(); i++)
        {
            final int idx = suppIndices.get(i);
            final SAMRecord supp = records.get(idx);
            final LiftBackResult suppRes = resolved[idx];
            if(suppRes == null || suppRes.finalCigar() == null
                    || suppRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
                continue;
            suppDtos.add(new RescueSupplementary(
                    i, suppRes.finalChrom(), !supp.getReadNegativeStrandFlag(),
                    suppRes.finalPos(), suppRes.finalCigar(), supp.getMappingQuality()));
        }

        final RescueCandidate candidate = new RescueCandidate(
                primaryRes.finalChrom(), !primary.getReadNegativeStrandFlag(), primary.getReadLength(),
                primaryRes.finalPos(), primaryRes.finalCigar(), primary.getMappingQuality(), suppDtos,
                primary.getReadBases(), mateHintIntrons);

        final RescueResult result = mJunctionRescueResolver.resolve(candidate);
        if(!result.merged())
            return null;

        // rewrite the primary's LiftBackResult with the merged (now-spliced) cigar + start; the start may
        // shift for a left-extend, and MAPQ caps at RESCUED_MAPQ_CAP to mark it a constructed alignment.
        // mark merged supps for drop.
        resolved[primaryIdx] = primaryRes.withRevisedCigar(
                result.MergedStart, result.MergedCigar, true,
                Math.min(primaryRes.updatedMapq(), RESCUED_MAPQ_CAP), "rescued-via-supp");

        if(introducedIntronsOut != null && result.IntroducedIntrons != null)
            introducedIntronsOut.addAll(result.IntroducedIntrons);

        final boolean[] dropped = new boolean[records.size()];
        for(Integer dtoIdx : result.DroppedSupplementaryIndices)
            dropped[suppIndices.get(dtoIdx)] = true;
        return dropped;
    }

    // Runs after applyJunctionRescue so rescue's lookups see bwa's original cigar; the extender then
    // cleans up tail-trim residual the rescue couldn't merge. Skipped on rescue-merged primaries.
    private void applyTerminalJunctionCollapse(
            final List<SAMRecord> records, final LiftBackResult[] resolved, final SAMRecord primary,
            final boolean[] droppedByRescue)
    {
        if(mTerminalJunctionCollapser == null || primary == null || primary.getReadUnmappedFlag())
            return;

        final int primaryIdx = indexOfPrimary(records, primary);
        if(primaryIdx < 0 || (droppedByRescue != null && droppedByRescue[primaryIdx]))
            return;

        final LiftBackResult primaryRes = resolved[primaryIdx];
        if(primaryRes == null || primaryRes.finalCigar() == null
                || primaryRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
            return;
        if(!primaryRes.hasNCigar())
            return;

        final TerminalCollapseResult collapse = mTerminalJunctionCollapser.tryCollapse(
                primaryRes.finalChrom(), primaryRes.finalPos(), primaryRes.finalCigar(), primary.getReadBases());

        if(!collapse.Collapsed)
            return;

        resolved[primaryIdx] = primaryRes.withRevisedCigar(
                collapse.NewStart, collapse.NewCigar, collapse.NewCigar.indexOf('N') >= 0,
                primaryRes.updatedMapq(), "terminal-junction-collapsed");
    }

    private void applyTailExtension(
            final List<SAMRecord> records, final LiftBackResult[] resolved, final SAMRecord primary,
            final boolean[] droppedByRescue)
    {
        if(mSoftclipTailExtender == null || primary == null || primary.getReadUnmappedFlag())
            return;

        final int primaryIdx = indexOfPrimary(records, primary);
        if(primaryIdx < 0)
            return;

        if(droppedByRescue != null && droppedByRescue[primaryIdx])
            return;

        final LiftBackResult primaryRes = resolved[primaryIdx];
        if(primaryRes == null || primaryRes.finalCigar() == null
                || primaryRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
            return;

        final TailExtensionResult extension = mSoftclipTailExtender.tryExtend(
                primaryRes.finalChrom(), primaryRes.finalPos(), primaryRes.finalCigar(),
                primary.getReadBases());

        if(!extension.Extended)
            return;

        resolved[primaryIdx] = primaryRes.withRevisedCigar(
                extension.NewStart, extension.NewCigar, primaryRes.hasNCigar(),
                primaryRes.updatedMapq(), "tail-extended");
    }

    private void applyJunctionCanonicalization(
            final List<SAMRecord> records, final LiftBackResult[] resolved, final SAMRecord primary,
            final boolean[] droppedByRescue)
    {
        if(mJunctionCanonicalizer == null || primary == null || primary.getReadUnmappedFlag())
            return;

        final int primaryIdx = indexOfPrimary(records, primary);
        if(primaryIdx < 0 || (droppedByRescue != null && droppedByRescue[primaryIdx]))
            return;

        final LiftBackResult primaryRes = resolved[primaryIdx];
        if(primaryRes == null || primaryRes.finalCigar() == null
                || primaryRes.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
            return;
        if(!primaryRes.hasNCigar())
            return;

        final JunctionCanonicalizationResult canon = mJunctionCanonicalizer.tryCanonicalize(
                primaryRes.finalChrom(), primaryRes.finalPos(), primaryRes.finalCigar(), primary.getReadBases());

        if(!canon.Changed)
            return;

        resolved[primaryIdx] = primaryRes.withRevisedCigar(
                primaryRes.finalPos(), canon.NewCigar, primaryRes.hasNCigar(),
                primaryRes.updatedMapq(), "junction-canonicalized");
    }

    // Dedup key for supplementaries: lifted (chrom, pos, cigar) + lifted strand. Opposite-strand
    // placements at the same coords/cigar are kept as distinct records.
    private static String dedupKey(final LiftBackResult result, final SAMRecord record)
    {
        return result.finalChrom() + ":" + result.finalPos() + ":" + result.finalCigar()
                + ":" + (result.negativeStrand() ? '-' : '+');
    }

    private void applyAndWriteRecord(
            final SAMRecord record, final LiftBackResult result, final int nh, final boolean primaryPostProcessed,
            final LiftedMateInfoCache liftedMateInfoCache, final EmitSink sink)
    {
        mStats.record(record, result);

        LiftBackRecordOps.applyResultToRecord(record, result, liftedMateInfoCache);

        final String rewrittenSa = rewriteSaTag(record.getStringAttribute(SA_ATTRIBUTE), mResolver);

        // A record that stays supplementary must carry an SA tag -- REDUX dedup (FragmentCoords.fromRead)
        // reads the primary's coords from it. If every SA entry failed to lift the supplementary is orphaned
        // in genomic space and cannot be represented, so drop it rather than emit a supp with a null SA.
        if(rewrittenSa == null && record.getSupplementaryAlignmentFlag())
        {
            mStats.recordOrphanSuppDropped();
            return;
        }

        record.setAttribute(SA_ATTRIBUTE, rewrittenSa);
        patchMateFields(record, liftedMateInfoCache);

        if(result.category() != LiftBackCategory.UNMAPPED)
            record.setAttribute("NH", nh);

        if(mUnmapAboveNh > 0 && nh > mUnmapAboveNh
                && !record.getReadUnmappedFlag()
                && !record.getSupplementaryAlignmentFlag())
        {
            markPrimaryUnmapped(record);
        }

        if(mUnmapBelowMapq > 0
                && !record.getReadUnmappedFlag()
                && !record.getSupplementaryAlignmentFlag()
                && record.getMappingQuality() < mUnmapBelowMapq)
        {
            markPrimaryUnmapped(record);
        }

        // residual short-anchor primary: bwa scored it below the default -T 30 floor and liftback couldn't
        // improve it (not rescued / extended / collapsed). Unmap it so the noise alignment doesn't pollute
        // the locus. Gated on rescue running, matching the supp drop's -T 19 rationale.
        if(mJunctionRescueResolver != null && !primaryPostProcessed
                && !record.getReadUnmappedFlag()
                && !record.getSupplementaryAlignmentFlag())
        {
            final Integer alignmentScore = record.getIntegerAttribute(AS_TAG);
            if(alignmentScore != null && alignmentScore < PRIMARY_AS_UNMAP_THRESHOLD)
            {
                markPrimaryUnmapped(record);
                mStats.recordLowAsPrimaryUnmapped();
            }
        }

        refreshNmDropMd(record, mRefSource);
        sink.emit(record, result);
    }
}
