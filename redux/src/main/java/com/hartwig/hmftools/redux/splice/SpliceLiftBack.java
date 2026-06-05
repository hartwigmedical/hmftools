package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.splice.MateFieldPatcher.patchMateFields;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.SA_ATTRIBUTE;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.rewriteSaTag;
import static com.hartwig.hmftools.redux.splice.SpliceCommon.ALT_CONTIG_SUFFIX;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.redux.splice.rescue.AnnotatedJunctionIndex;
import com.hartwig.hmftools.redux.splice.rescue.AnnotatedJunctionLoader;
import com.hartwig.hmftools.redux.splice.rescue.ChrIntron;
import com.hartwig.hmftools.redux.splice.rescue.JunctionRescueResolver;
import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;
import com.hartwig.hmftools.redux.splice.rescue.RescueCandidate;
import com.hartwig.hmftools.redux.splice.rescue.RescueConfig;
import com.hartwig.hmftools.redux.splice.rescue.RescueResult;
import com.hartwig.hmftools.redux.splice.rescue.RescueStatistics;
import com.hartwig.hmftools.redux.splice.rescue.RescueSupplementary;
import com.hartwig.hmftools.redux.splice.tailextend.SoftclipTailExtender;
import com.hartwig.hmftools.redux.splice.tailextend.TerminalCollapseResult;
import com.hartwig.hmftools.redux.splice.tailextend.TerminalMicroJunctionCollapser;
import com.hartwig.hmftools.redux.splice.tailextend.TailExtensionConfig;
import com.hartwig.hmftools.redux.splice.tailextend.TailExtensionResult;
import com.hartwig.hmftools.redux.splice.tailextend.TailExtensionStatistics;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;

// Two passes over the bwa-mem2 BAM aligned against ref + tx contigs:
//   pass 1: build a LiftedMateInfoCache so each record's partner-mate fields can be rewritten in genomic coords.
//   pass 2: resolve + apply + write a coordinate-sorted BAM with mate fields patched.
// 1:1 record contract: every input record produces exactly one output record + one TSV-A row + N TSV-B rows.
public class SpliceLiftBack
{
    private static final String XA_TAG = "XA";
    private static final String AS_TAG = "AS";

    // bwa-mem2's default -T (minimum alignment score) is 30. The RNA splice run uses -T 19 deliberately
    // to surface short-anchor supplementary records that JunctionRescueResolver can merge across annotated
    // junctions. Supps in the [19, 30) AS band that survive the rescue pass are residual noise -- they
    // were emitted only because of the lowered -T and were not promoted to spliced primaries. Drop them
    // from the output BAM after rescue + tail-extend have had their chance.
    static final int SUPP_AS_DROP_THRESHOLD = 30;

    private final SpliceLiftBackConfig mConfig;

    // name-sorted intermediate so primary + supps are contiguous in the stream.
    private String mWorkingInputBam;

    // optional rescue resolver: merges primary + supp across annotated junctions when bwa emits the
    // spliced fragment as two records. Null when -rescue_via_supp is not set.
    private JunctionRescueResolver mJunctionRescueResolver;

    // null when -extend_softclip_tails is not set.
    private SoftclipTailExtender mSoftclipTailExtender;
    private TerminalMicroJunctionCollapser mTerminalJunctionCollapser;

    // Terminal anchor at or below this length is a candidate for the ref-parsimony collapse: the
    // tx-contig walk fabricates a junction with a 1-2bp exon anchor at a read end when the read runs
    // a base or two past an exon boundary on the concatenated contig.
    private static final int TERMINAL_MICRO_JUNCTION_MAX_ANCHOR = 2;

    public SpliceLiftBack(final ConfigBuilder configBuilder)
    {
        mConfig = new SpliceLiftBackConfig(configBuilder);
    }

    public void run()
    {
        final long startTimeMs = System.currentTimeMillis();

        mWorkingInputBam = nameSortInput();

        final List<ContigEntry> contigEntries = ContigSidecar.read(mConfig.ContigSidecarFile);
        RD_LOGGER.info("loaded {} contig entries", contigEntries.size());

        ExonRegionIndex exonIndex = null;
        if(mConfig.EnsemblDataDir != null)
        {
            try
            {
                exonIndex = ExonRegionIndex.load(mConfig.EnsemblDataDir);
                RD_LOGGER.info("loaded annotated-exon index from {}", mConfig.EnsemblDataDir);
            }
            catch(IOException e)
            {
                throw new IllegalStateException("failed to load annotated exons: " + e, e);
            }
        }

        final LiftBackResolver resolver = new LiftBackResolver(contigEntries, exonIndex);
        validateBamAgainstSidecar(resolver.contigNames());

        if(mConfig.RescueViaSupp || mConfig.ExtendSoftclipTails)
        {
            try
            {
                final boolean needAnnotation = mConfig.RescueViaSupp || mConfig.ExtendSoftclipTails;
                final AnnotatedJunctionIndex index;
                if(needAnnotation)
                {
                    final Set<ChrIntron> annotated = AnnotatedJunctionLoader.load(mConfig.EnsemblDataDir);
                    index = new AnnotatedJunctionIndex(annotated);
                    RD_LOGGER.info("loaded {} annotated junctions from {}",
                            annotated.size(), mConfig.EnsemblDataDir);
                }
                else
                {
                    index = null;
                }

                final RefSequenceSource refSource = openRefSource(mConfig.RefGenomeFile);

                if(mConfig.RescueViaSupp)
                {
                    mJunctionRescueResolver = new JunctionRescueResolver(
                            index, refSource, RescueConfig.enabledDefaults());
                    RD_LOGGER.info("rescue-via-supp enabled (ref-verify: {})",
                            refSource != null ? "enabled" : "disabled");
                }

                if(mConfig.ExtendSoftclipTails)
                {
                    mSoftclipTailExtender = new SoftclipTailExtender(
                            refSource, index, TailExtensionConfig.enabledDefaults());
                    RD_LOGGER.info("extend-softclip-tails enabled (junction guard: {})",
                            index != null ? "enabled" : "disabled");
                }

                // Runs whenever ref is available: collapses spurious tx-contig terminal micro-junctions
                // (a 1-2bp exon anchor across an intron at a read end) when the bases align contiguously.
                if(refSource != null)
                {
                    mTerminalJunctionCollapser = new TerminalMicroJunctionCollapser(
                            refSource, TERMINAL_MICRO_JUNCTION_MAX_ANCHOR);
                    RD_LOGGER.info("terminal micro-junction collapse enabled (max anchor: {})",
                            TERMINAL_MICRO_JUNCTION_MAX_ANCHOR);
                }
            }
            catch(IOException e)
            {
                throw new IllegalStateException("failed to load annotated junctions: " + e, e);
            }
        }

        final LiftedMateInfoCache liftedMateInfoCache = buildLiftedMateInfoCache(resolver);

        emitLiftedBam(resolver, liftedMateInfoCache);

        if(mJunctionRescueResolver != null)
            logRescueStats(mJunctionRescueResolver.statistics());

        if(mSoftclipTailExtender != null)
            logTailExtensionStats(mSoftclipTailExtender.statistics());

        if(mTerminalJunctionCollapser != null)
            RD_LOGGER.info("terminal micro-junction collapse: leading={} trailing={}",
                    mTerminalJunctionCollapser.collapsedLeading(), mTerminalJunctionCollapser.collapsedTrailing());

        try
        {
            Files.deleteIfExists(Paths.get(mWorkingInputBam));
        }
        catch(IOException e)
        {
            RD_LOGGER.warn("could not delete name-sorted intermediate {}: {}", mWorkingInputBam, e.toString());
        }

        RD_LOGGER.info("SpliceLiftBack complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private String nameSortInput()
    {
        if(mConfig.BamToolPath == null)
            throw new IllegalStateException("name-sort pre-pass requires -bamtool");

        final String nameSortedBam = mConfig.OutputDir + "splice_lift_back.name_sorted_input.bam";
        RD_LOGGER.info("name-sorting input via {}: {}", fromPath(mConfig.BamToolPath), nameSortedBam);

        final long startTimeMs = System.currentTimeMillis();
        final java.util.List<String> command = new ArrayList<>();
        command.add(mConfig.BamToolPath);
        command.add("sort");
        command.add("-n");
        if(mConfig.Threads > 1)
        {
            command.add("-@");
            command.add(String.valueOf(mConfig.Threads));
        }
        if(fromPath(mConfig.BamToolPath) == BamToolName.SAMTOOLS)
        {
            command.add("-O");
            command.add("bam");
        }
        command.add(mConfig.InputBam);
        command.add("-o");
        command.add(nameSortedBam);

        try
        {
            final Process process = new ProcessBuilder(command).redirectErrorStream(true).start();
            final int exitCode = process.waitFor();
            if(exitCode != 0)
                throw new RuntimeException("name-sort exit code " + exitCode);
        }
        catch(InterruptedException | IOException e)
        {
            throw new RuntimeException("name-sort failed: " + e, e);
        }

        RD_LOGGER.info("name-sort complete, mins({})", runTimeMinsStr(startTimeMs));
        return nameSortedBam;
    }

    // fails fast on a BAM/sidecar mismatch — root cause of a prior silent failure where alt contigs
    // missing from the sidecar fell through as ref pass-through and leaked _tx names into mate fields.
    private void validateBamAgainstSidecar(final Set<String> sidecarContigs)
    {
        Set<String> missing = new TreeSet<>();
        try(SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mWorkingInputBam)))
        {
            for(SAMSequenceRecord sq : reader.getFileHeader().getSequenceDictionary().getSequences())
            {
                String name = sq.getSequenceName();
                if(name.endsWith(ALT_CONTIG_SUFFIX) && !sidecarContigs.contains(name))
                    missing.add(name);
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read BAM header: " + mWorkingInputBam, e);
        }

        if(!missing.isEmpty())
        {
            throw new IllegalStateException(String.format(
                    "BAM/sidecar mismatch: %d alt contig(s) in BAM @SQ are absent from sidecar %s — first few: %s",
                    missing.size(), mConfig.ContigSidecarFile,
                    missing.stream().limit(5).collect(Collectors.joining(","))));
        }
    }

    private LiftedMateInfoCache buildLiftedMateInfoCache(final LiftBackResolver resolver)
    {
        final long startTimeMs = System.currentTimeMillis();
        final LiftedMateInfoCache cache = new LiftedMateInfoCache();
        int primaryRecordsCached = 0;

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mWorkingInputBam)))
        {
            final SAMRecordIterator iter = samReader.iterator();
            while(iter.hasNext())
            {
                final SAMRecord record = iter.next();

                if(!record.getReadPairedFlag())
                    continue;
                if(record.isSecondaryOrSupplementary())
                    continue;

                final LiftBackResult result = resolver.resolve(record);
                final LiftedMateInfo info = toLiftedMateInfo(record, result, mConfig.UnmapAboveNh, mConfig.UnmapBelowMapq);
                cache.recordPrimaryAlignment(record.getReadName(), record.getFirstOfPairFlag(), info);
                ++primaryRecordsCached;
            }
        }
        catch(IOException e)
        {
            RD_LOGGER.error("pass 1 BAM read failed", e);
            throw new RuntimeException(e);
        }

        RD_LOGGER.info("pass 1 cached {} primary alignments across {} read pairs, mins({})",
                primaryRecordsCached, cache.size(), runTimeMinsStr(startTimeMs));
        return cache;
    }

    private void emitLiftedBam(final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache)
    {
        final long startTimeMs = System.currentTimeMillis();
        final LiftBackStats stats = new LiftBackStats();

        final String unsortedBam = mConfig.formUnsortedBam();

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mWorkingInputBam));
            LiftBackWriter writer = new LiftBackWriter(mConfig.formTsvAFile(), mConfig.formTsvBFile()))
        {
            final SAMFileHeader header = samReader.getFileHeader().clone();
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            stripAltContigsFromHeader(header);

            try(SAMFileWriter samWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, new File(unsortedBam)))
            {
                final SAMRecordIterator iter = samReader.iterator();
                streamGrouped(iter, samWriter, writer, resolver, liftedMateInfoCache, stats);
            }

            stats.writeSummary(mConfig.formSummaryFile());
        }
        catch(IOException e)
        {
            RD_LOGGER.error("pass 2 BAM/TSV I/O failed", e);
            throw new RuntimeException(e);
        }

        stats.logSummary();
        RD_LOGGER.info("pass 2 emit complete, mins({})", runTimeMinsStr(startTimeMs));

        sortAndIndex(unsortedBam, mConfig.formOutputBam());
    }

    // Buffers records sharing a read name so primary + supps process as one group. This lets rescue
    // see all split-read components and lets /2 hint /1's introns (and vice versa). Supplementaries
    // are split-read components, not alts — alts arrive on the primary's XA tag and are handled by
    // the resolver per-record.
    private void streamGrouped(
            final SAMRecordIterator iter, final SAMFileWriter samWriter, final LiftBackWriter writer,
            final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache, final LiftBackStats stats)
    {
        final List<SAMRecord> buffer = new ArrayList<>();
        String currentName = null;
        while(iter.hasNext())
        {
            final SAMRecord record = iter.next();
            final String name = record.getReadName();
            if(currentName != null && !name.equals(currentName))
            {
                flushReadNameGroup(buffer, samWriter, writer, resolver, liftedMateInfoCache, stats);
                buffer.clear();
            }
            buffer.add(record);
            currentName = name;
        }
        if(!buffer.isEmpty())
        {
            flushReadNameGroup(buffer, samWriter, writer, resolver, liftedMateInfoCache, stats);
        }
    }

    private void flushReadNameGroup(
            final List<SAMRecord> group, final SAMFileWriter samWriter, final LiftBackWriter writer,
            final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache, final LiftBackStats stats)
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
        final int unmapAboveNh = mConfig.UnmapAboveNh;
        final int unmapBelowMapq = mConfig.UnmapBelowMapq;
        final MateDecision m1 = decideMateGroup(firstOfPair, resolver, java.util.Collections.emptyList());
        refreshMateInfoCache(firstOfPair, m1, liftedMateInfoCache, unmapAboveNh, unmapBelowMapq);
        final MateDecision m2 = decideMateGroup(secondOfPair, resolver, m1.IntroducedIntrons);
        refreshMateInfoCache(secondOfPair, m2, liftedMateInfoCache, unmapAboveNh, unmapBelowMapq);
        final MateDecision m1Final = (m1.IntroducedIntrons.isEmpty() && !m2.IntroducedIntrons.isEmpty())
                ? decideMateGroup(firstOfPair, resolver, m2.IntroducedIntrons)
                : m1;
        if(m1Final != m1)
            refreshMateInfoCache(firstOfPair, m1Final, liftedMateInfoCache, unmapAboveNh, unmapBelowMapq);

        writeMateGroup(firstOfPair, m1Final, samWriter, writer, resolver, liftedMateInfoCache, stats);
        writeMateGroup(secondOfPair, m2, samWriter, writer, resolver, liftedMateInfoCache, stats);
    }

    // Push the post-rescue primary back into the mate cache so the partner mate's MC tag uses
    // the merged/extended cigar instead of the pre-rescue one pass 1 seeded.
    private static void refreshMateInfoCache(
            final List<SAMRecord> mateRecords, final MateDecision decision, final LiftedMateInfoCache cache,
            final int unmapAboveNh, final int unmapBelowMapq)
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
        final LiftedMateInfo refreshed = toLiftedMateInfo(
                primary, decision.PrimaryResult, unmapAboveNh, unmapBelowMapq);
        cache.recordPrimaryAlignment(primary.getReadName(), primary.getFirstOfPairFlag(), refreshed);
    }

    // Per-mate decision: resolves every record, runs the rescue resolver, returns the lifted
    // results + per-record drop flags + introns introduced by rescue (for mate hinting). Does not
    // write to BAM — that's the second phase (writeMateGroup).
    private static final class MateDecision
    {
        final LiftBackResult[] Resolved;
        final boolean[] DroppedByRescue;
        final LiftBackResult PrimaryResult;
        final List<ChrIntron> IntroducedIntrons;

        MateDecision(final LiftBackResult[] resolved, final boolean[] droppedByRescue,
                final LiftBackResult primaryResult,
                final List<ChrIntron> introducedIntrons)
        {
            Resolved = resolved;
            DroppedByRescue = droppedByRescue;
            PrimaryResult = primaryResult;
            IntroducedIntrons = introducedIntrons != null ? introducedIntrons : java.util.Collections.emptyList();
        }

        static MateDecision empty()
        {
            return new MateDecision(new LiftBackResult[0], null, null, java.util.Collections.emptyList());
        }
    }

    private MateDecision decideMateGroup(
            final List<SAMRecord> records, final LiftBackResolver resolver,
            final List<ChrIntron> mateHintIntrons)
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

        final LiftBackResult primaryResult = primary != null ? resolver.resolve(primary) : null;

        final LiftBackResult[] resolved = new LiftBackResult[records.size()];
        for(int i = 0; i < records.size(); i++)
        {
            final SAMRecord record = records.get(i);
            resolved[i] = (record == primary && primaryResult != null) ? primaryResult : resolver.resolve(record);
        }

        final List<ChrIntron> introduced = new ArrayList<>();
        final boolean[] droppedByRescue = applyJunctionRescue(records, resolved, primary, mateHintIntrons, introduced);

        // Collapse a spurious tx-contig terminal micro-junction before tail extension, so the extender
        // sees the corrected (contiguous) cigar rather than the fabricated junction.
        applyTerminalJunctionCollapse(records, resolved, primary, droppedByRescue);

        applyTailExtension(records, resolved, primary, droppedByRescue);

        return new MateDecision(resolved, droppedByRescue,
                resolved.length > 0 && primary != null ? resolved[indexOfPrimary(records, primary)] : primaryResult,
                introduced);
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
            final SAMFileWriter samWriter, final LiftBackWriter writer,
            final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache,
            final LiftBackStats stats)
    {
        if(records.isEmpty())
            return;

        final LiftBackResult[] resolved = decision.Resolved;
        final boolean[] droppedByRescue = decision.DroppedByRescue;
        SAMRecord primary = null;
        for(SAMRecord r : records)
            if(!r.getSupplementaryAlignmentFlag()) { primary = r; break; }

        // Dedup supplementaries that lift to the same (chrom, pos, cigar, strand) — bwa can emit the
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
                    stats.recordLowAsSuppDropped();
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
                stats.record(record, result);
                continue;
            }
            applyAndWriteRecord(record, result, nh, samWriter, writer, resolver, liftedMateInfoCache, stats);
        }
    }

    // builds a RescueCandidate from the post-lift primary + its supplementary records and applies
    // the merge if accepted. Returns a parallel array marking which records were absorbed by the
    // merge and should be dropped from emission. Returns null when rescue is disabled / no-op.
    // introducedIntronsOut is appended-to (caller-provided list) so the partner mate can use them
    // as hints in the second pass.
    private boolean[] applyJunctionRescue(
            final List<SAMRecord> records, final LiftBackResult[] resolved, final SAMRecord primary,
            final List<ChrIntron> mateHintIntrons, final List<ChrIntron> introducedIntronsOut)
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

        // rewrite the primary's LiftBackResult with the merged cigar + start; mark merged supps for drop.
        resolved[primaryIdx] = mergedPrimaryResult(primaryRes, result.MergedCigar, result.MergedStart);

        if(introducedIntronsOut != null && result.IntroducedIntrons != null)
            introducedIntronsOut.addAll(result.IntroducedIntrons);

        final boolean[] dropped = new boolean[records.size()];
        for(Integer dtoIdx : result.DroppedSupplementaryIndices)
            dropped[suppIndices.get(dtoIdx)] = true;
        return dropped;
    }

    // Runs after applyJunctionRescue so rescue's lookups see bwa's original cigar; the extender
    // then cleans up tail-trim residual the rescue couldn't merge. Skipped on rescue-merged
    // primaries — their terminal softclips are already gone.
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

        resolved[primaryIdx] = terminalCollapsedResult(primaryRes, collapse);
    }

    private static LiftBackResult terminalCollapsedResult(
            final LiftBackResult original, final TerminalCollapseResult collapse)
    {
        final boolean stillHasN = collapse.NewCigar.indexOf('N') >= 0;
        return new LiftBackResult(
                original.category(), original.comp(), original.role(),
                original.finalChrom(), collapse.NewStart, collapse.NewCigar,
                original.negativeStrand(), stillHasN,
                original.inputMapq(), original.updatedMapq(),
                original.numXaAlts(), original.numRefAlts(), original.numTxAlts(),
                original.numLoci(), original.numDistinctCigarsAtPrimaryLocus(),
                original.txHasNCigar(), original.txSoftClipAtBoundary(),
                original.refSoftClipped(), original.refFullMatch(),
                original.geneIds(),
                appendNote(original.notes(), "terminal-junction-collapsed"),
                original.transcriptStrand(),
                original.liftedAlignments());
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

        resolved[primaryIdx] = tailExtendedResult(primaryRes, extension);
    }

    private static LiftBackResult tailExtendedResult(
            final LiftBackResult original, final TailExtensionResult extension)
    {
        return new LiftBackResult(
                original.category(), original.comp(), original.role(),
                original.finalChrom(), extension.NewStart, extension.NewCigar,
                original.negativeStrand(), original.hasNCigar(),
                original.inputMapq(), original.updatedMapq(),
                original.numXaAlts(), original.numRefAlts(), original.numTxAlts(),
                original.numLoci(), original.numDistinctCigarsAtPrimaryLocus(),
                original.txHasNCigar(), original.txSoftClipAtBoundary(),
                original.refSoftClipped(), original.refFullMatch(),
                original.geneIds(),
                appendNote(original.notes(), "tail-extended"),
                original.transcriptStrand(),
                original.liftedAlignments());
    }

    // BWA emits MAPQ=60 for a clean unique alignment. Rescued primaries are constructed by us, not
    // directly emitted by BWA, so we cap at 55 to give downstream tools a signal that this record
    // came from a primary+supp merge rather than a single confident BWA placement. If the original
    // MAPQ was already below 55 (e.g. MAPQ=0 multi-mapper), keep it — capping never goes UP.
    private static final int RESCUED_MAPQ_CAP = 55;

    // constructs an updated LiftBackResult preserving every field except finalPos/finalCigar/
    // hasNCigar/updatedMapq. The rescue changes the cigar shape (now spliced), possibly the start
    // position (for left-extends), and caps MAPQ to RESCUED_MAPQ_CAP to mark it as a constructed
    // alignment.
    private static LiftBackResult mergedPrimaryResult(
            final LiftBackResult original, final String mergedCigar, final int mergedStart)
    {
        final int rescuedMapq = Math.min(original.updatedMapq(), RESCUED_MAPQ_CAP);
        return new LiftBackResult(
                original.category(), original.comp(), original.role(),
                original.finalChrom(), mergedStart, mergedCigar,
                original.negativeStrand(), true,
                original.inputMapq(), rescuedMapq,
                original.numXaAlts(), original.numRefAlts(), original.numTxAlts(),
                original.numLoci(), original.numDistinctCigarsAtPrimaryLocus(),
                original.txHasNCigar(), original.txSoftClipAtBoundary(),
                original.refSoftClipped(), original.refFullMatch(),
                original.geneIds(),
                appendNote(original.notes(), "rescued-via-supp"),
                original.transcriptStrand(),
                original.liftedAlignments());
    }

    private static String appendNote(final String existing, final String note)
    {
        if(existing == null || existing.isEmpty())
            return note;
        return existing + ";" + note;
    }

    // Opens the reference FASTA for ref-verify lookups. Returns a RefSequenceSource that wraps
    // htsjdk's IndexedFastaSequenceFile. Synchronized internally because the indexed reader is
    // not thread-safe and the rescue path may be called concurrently from per-partition workers.
    private static RefSequenceSource openRefSource(final String refGenomeFile)
    {
        if(refGenomeFile == null)
            return null;
        try
        {
            final IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(new File(refGenomeFile));
            return (chromosome, posStart, posEnd) ->
            {
                synchronized(fasta)
                {
                    try
                    {
                        final ReferenceSequence seq = fasta.getSubsequenceAt(chromosome, posStart, posEnd);
                        return seq != null ? seq.getBases() : null;
                    }
                    catch(Exception e)
                    {
                        return null;
                    }
                }
            };
        }
        catch(Exception e)
        {
            RD_LOGGER.warn("could not open ref FASTA {} for ref-verify: {}", refGenomeFile, e.toString());
            return null;
        }
    }

    private static void logTailExtensionStats(final TailExtensionStatistics stats)
    {
        RD_LOGGER.info("extend-softclip-tails summary: evaluated={} extended={} basesLead={} basesTrail={}",
                stats.recordsEvaluated(), stats.recordsExtended(),
                stats.basesExtendedLead(), stats.basesExtendedTrail());
        if(stats.skippedNoRef() > 0)
            RD_LOGGER.info("extend-softclip-tails skipped no-ref: {}", stats.skippedNoRef());
        if(stats.skippedForJunctionGuard() > 0)
            RD_LOGGER.info("extend-softclip-tails skipped junction-guard: {}", stats.skippedForJunctionGuard());
        if(stats.skippedComplexShape() > 0)
            RD_LOGGER.info("extend-softclip-tails skipped complex-shape: {}", stats.skippedComplexShape());
        if(stats.rejectedTooManyMismatches() > 0)
            RD_LOGGER.info("extend-softclip-tails rejected too-many-mismatches: {}",
                    stats.rejectedTooManyMismatches());
    }

    private static void logRescueStats(final RescueStatistics stats)
    {
        RD_LOGGER.info("rescue-via-supp summary: candidates={} merged={} suppClamped={} (depth1={} depth2={} depth3={} depth4={})",
                stats.candidatesEvaluated(), stats.mergedTotal(), stats.suppClampApplied(),
                stats.mergedAtChainDepth(1), stats.mergedAtChainDepth(2),
                stats.mergedAtChainDepth(3), stats.mergedAtChainDepth(4));
        for(com.hartwig.hmftools.redux.splice.rescue.RescueRejectReason reason
                : com.hartwig.hmftools.redux.splice.rescue.RescueRejectReason.values())
        {
            final int count = stats.rejectCount(reason);
            if(count > 0)
                RD_LOGGER.info("rescue-via-supp reject {}: {}", reason, count);
        }
    }

    // Dedup key for supplementaries: lifted (chrom, pos, cigar) + lifted strand. Opposite-strand
    // placements at the same coords/cigar are kept as distinct records.
    private static String dedupKey(final LiftBackResult result, final SAMRecord record)
    {
        return result.finalChrom() + ":" + result.finalPos() + ":" + result.finalCigar()
                + ":" + (result.negativeStrand() ? '-' : '+');
    }

    private void applyAndWriteRecord(
            final SAMRecord record, final LiftBackResult result, final int nh, final SAMFileWriter samWriter,
            final LiftBackWriter writer, final LiftBackResolver resolver,
            final LiftedMateInfoCache liftedMateInfoCache, final LiftBackStats stats)
    {
        stats.record(record, result);

        applyResultToRecord(record, result, liftedMateInfoCache);
        record.setAttribute(SA_ATTRIBUTE, rewriteSaTag(record.getStringAttribute(SA_ATTRIBUTE), resolver));
        patchMateFields(record, liftedMateInfoCache);

        if(result.category() != LiftBackCategory.UNMAPPED)
            record.setAttribute("NH", nh);

        if(mConfig.UnmapAboveNh > 0 && nh > mConfig.UnmapAboveNh
                && !record.getReadUnmappedFlag()
                && !record.getSupplementaryAlignmentFlag())
        {
            markPrimaryUnmapped(record);
        }

        if(mConfig.UnmapBelowMapq > 0
                && !record.getReadUnmappedFlag()
                && !record.getSupplementaryAlignmentFlag()
                && record.getMappingQuality() < mConfig.UnmapBelowMapq)
        {
            markPrimaryUnmapped(record);
        }

        stripMdNm(record);
        samWriter.addAlignment(record);

        try
        {
            writer.write(record, result);
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    // MD/NM from the original tx-contig alignment are stale post-lift; downstream tools recompute them.
    private static void stripMdNm(final SAMRecord record)
    {
        record.setAttribute("MD", null);
        record.setAttribute("NM", null);
    }

    private void sortAndIndex(final String unsortedBam, final String sortedBam)
    {
        if(mConfig.BamToolPath == null)
        {
            RD_LOGGER.info("no -{} configured; leaving unsorted BAM at {}", BamToolName.BAMTOOL_PATH, unsortedBam);
            return;
        }

        final BamToolName toolName = fromPath(mConfig.BamToolPath);
        RD_LOGGER.info("sorting BAM via {}: {}", toolName, sortedBam);

        if(!BamOperations.sortBam(toolName, mConfig.BamToolPath, unsortedBam, sortedBam, mConfig.Threads))
            throw new RuntimeException("failed to sort BAM: " + unsortedBam);

        if(toolName == BamToolName.SAMTOOLS)
        {
            if(!BamOperations.indexBam(toolName, mConfig.BamToolPath, sortedBam, mConfig.Threads))
                throw new RuntimeException("failed to index BAM: " + sortedBam);
        }

        try
        {
            Files.deleteIfExists(Paths.get(unsortedBam));
        }
        catch(IOException e)
        {
            RD_LOGGER.warn("could not delete intermediate {}: {}", unsortedBam, e.toString());
        }
    }

    private static void stripAltContigsFromHeader(final SAMFileHeader header)
    {
        final SAMSequenceDictionary dict = header.getSequenceDictionary();
        final List<SAMSequenceRecord> kept = new ArrayList<>();
        int dropped = 0;

        for(final SAMSequenceRecord seq : dict.getSequences())
        {
            if(seq.getSequenceName().endsWith(ALT_CONTIG_SUFFIX))
                ++dropped;
            else
                kept.add(seq);
        }

        header.setSequenceDictionary(new SAMSequenceDictionary(kept));
        RD_LOGGER.info("dropped {} alt contig @SQ entries from output header", dropped);
    }

    // Supplementaries whose own lift failed are mirrored onto their primary's lifted coords (keeping
    // the 0x800 flag) rather than marked unmapped — htsjdk's validator rejects 0x4+0x800.
    static void applyResultToRecord(
            final SAMRecord record, final LiftBackResult result, final LiftedMateInfoCache liftedMateInfoCache)
    {
        switch(result.category())
        {
            case UNMAPPED:
                return;

            case LIFT_FAILED:
                if(record.isSecondaryOrSupplementary() && mirrorOwnPrimaryOntoFailedSupp(record, liftedMateInfoCache))
                    return;

                record.setReadUnmappedFlag(true);
                record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
                record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
                record.setMappingQuality(0);
                record.setAttribute(XA_TAG, null);
                record.setAttribute(SA_ATTRIBUTE, null);
                return;

            default:
                final Cigar liftedCigar = TextCigarCodec.decode(result.finalCigar());
                record.setReferenceName(result.finalChrom());
                record.setAlignmentStart(result.finalPos());
                record.setCigar(liftedCigar);
                record.setReadNegativeStrandFlag(result.negativeStrand());
                record.setMappingQuality(result.updatedMapq());
                record.setAttribute(XA_TAG, buildLiftedXa(result));
                // STAR-style XS:A:+/- on spliced records: downstream RNA tools (Isofox) rely on
                // transcript strand for stranded junction interpretation. bwa-mem2 emits XS:i: as
                // the sub-optimal alignment score on the same tag name — clear first so the SAM
                // tag-type bookkeeping isn't ambiguous. Only set XS:A when the lifted cigar has
                // an N AND we know the transcript strand (i.e. the primary came off a tx contig).
                // Ref-only N-cigars from rescue/tail-extend don't have a strand threaded yet —
                // they ship without XS rather than risk a wrong call.
                record.setAttribute("XS", null);
                if(result.hasNCigar() && result.transcriptStrand() != 0)
                    record.setAttribute("XS", result.transcriptStrand() > 0 ? Character.valueOf('+') : Character.valueOf('-'));
        }
    }

    private static boolean mirrorOwnPrimaryOntoFailedSupp(
            final SAMRecord record, final LiftedMateInfoCache liftedMateInfoCache)
    {
        if(!record.getReadPairedFlag())
            return false;

        final LiftedMateInfo ownPrimary = liftedMateInfoCache.getOwnPrimary(record.getReadName(), record.getFirstOfPairFlag());
        if(ownPrimary == null || ownPrimary.unmapped())
            return false;

        record.setReferenceName(ownPrimary.chromosome());
        record.setAlignmentStart(ownPrimary.alignmentStart());
        record.setCigarString(ownPrimary.liftedCigar());
        record.setMappingQuality(0);
        record.setAttribute(XA_TAG, null);
        // SA stays — rewriteSaTag translates tx-contig entries downstream and FragmentCoords requires it on supps
        return true;
    }

    static LiftedMateInfo toLiftedMateInfo(final SAMRecord record, final LiftBackResult result)
    {
        return toLiftedMateInfo(record, result, 0, 0);
    }

    // unmapAboveNh / unmapBelowMapq are forwarded so cache entries reflect threshold-unmap decisions
    // that applyAndWriteRecord will make later. Without this, a partner record's patched mate fields
    // would still claim "mate mapped" for a primary about to be threshold-unmapped.
    static LiftedMateInfo toLiftedMateInfo(
            final SAMRecord record, final LiftBackResult result, final int unmapAboveNh, final int unmapBelowMapq)
    {
        if(willBeUnmapped(result, unmapAboveNh, unmapBelowMapq))
            return LiftedMateInfo.UNMAPPED;

        final Cigar liftedCigar = TextCigarCodec.decode(result.finalCigar());
        final int alignmentEnd = result.finalPos() + liftedCigar.getReferenceLength() - 1;
        return LiftedMateInfo.mapped(result.finalChrom(), result.finalPos(), alignmentEnd, result.finalCigar(), result.negativeStrand());
    }

    // Mark a previously-mapped primary as unmapped, clearing every per-record tag the SAM spec invariant
    // would otherwise leave inconsistent. SA/XA/NH/MC all reference now-stale coords; ProperPair and
    // TLEN are meaningless once one end is unmapped. Package-private for direct unit testing.
    static void markPrimaryUnmapped(final SAMRecord record)
    {
        record.setReadUnmappedFlag(true);
        record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
        record.setMappingQuality(0);
        record.setAttribute(SA_ATTRIBUTE, null);
        record.setAttribute(XA_TAG, null);
        record.setAttribute("NH", null);
        record.setAttribute("MC", null);
        if(record.getReadPairedFlag())
        {
            record.setProperPairFlag(false);
            record.setInferredInsertSize(0);
        }
    }

    // Single source of truth for whether a primary's final state ends up unmapped, used both at cache-
    // build time (so partner records see the correct MateUnmapped) and at apply time (where we mutate
    // the record itself). Mirrors the threshold logic in applyAndWriteRecord.
    static boolean willBeUnmapped(final LiftBackResult result, final int unmapAboveNh, final int unmapBelowMapq)
    {
        if(result.category() == LiftBackCategory.UNMAPPED || result.category() == LiftBackCategory.LIFT_FAILED)
            return true;
        final int nh = Math.max(result.numLoci(), 1);
        if(unmapAboveNh > 0 && nh > unmapAboveNh)
            return true;
        if(unmapBelowMapq > 0 && result.updatedMapq() < unmapBelowMapq)
            return true;
        return false;
    }

    static String buildLiftedXa(final LiftBackResult result)
    {
        final String xa = result.liftedAlignments().stream()
                .filter(la -> !la.IsPrimaryChoice && !la.Dropped)
                .map(SpliceLiftBack::formatXaEntry)
                .collect(Collectors.joining());
        return xa.isEmpty() ? null : xa;
    }

    private static String formatXaEntry(final LiftedAlignment la)
    {
        return la.LiftedChrom + ',' + (la.ForwardStrand ? '+' : '-') + la.LiftedPos + ','
                + la.LiftedCigar + ',' + la.NumMismatches + ';';
    }

    public static void main(final String[] args)
    {
        final ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SpliceLiftBackConfig.addConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        new SpliceLiftBack(configBuilder).run();
    }
}
