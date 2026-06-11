package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.splice.SpliceCommon.ALT_CONTIG_SUFFIX;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
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
import com.hartwig.hmftools.redux.splice.rescue.RescueConfig;
import com.hartwig.hmftools.redux.splice.rescue.RescueStatistics;
import com.hartwig.hmftools.redux.splice.tailextend.SoftclipTailExtender;
import com.hartwig.hmftools.redux.splice.tailextend.TerminalMicroJunctionCollapser;
import com.hartwig.hmftools.redux.splice.tailextend.TailExtensionConfig;
import com.hartwig.hmftools.redux.splice.tailextend.TailExtensionStatistics;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// Standalone driver over the bwa-mem2 BAM aligned against ref + tx contigs.
// Pass 1 builds LiftedMateInfoCache; pass 2 streams name-grouped records through LiftBackGroupProcessor.
// Per-record logic is shared with SpliceLiftBackStage (concurrent REDUX path).
public class SpliceLiftBack
{
    private final SpliceLiftBackConfig mConfig;

    // name-sorted intermediate so primary + supps are contiguous in the stream.
    private String mWorkingInputBam;

    // null when -rescue_via_supp is not set.
    private JunctionRescueResolver mJunctionRescueResolver;

    // null when -extend_softclip_tails is not set.
    private SoftclipTailExtender mSoftclipTailExtender;
    private TerminalMicroJunctionCollapser mTerminalJunctionCollapser;
    private JunctionCanonicalizer mJunctionCanonicalizer;

    // null when neither rescue nor extend is enabled.
    private RefSequenceSource mRefSource;

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
                mRefSource = refSource;

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

                if(refSource != null)
                {
                    mTerminalJunctionCollapser = new TerminalMicroJunctionCollapser(
                            refSource, SpliceCommon.MIN_JUNCTION_ANCHOR);
                    RD_LOGGER.info("terminal micro-junction collapse enabled (min junction anchor: {})",
                            SpliceCommon.MIN_JUNCTION_ANCHOR);

                    mJunctionCanonicalizer = new JunctionCanonicalizer(
                            refSource, JunctionCanonicalizer.DEFAULT_MAX_SHIFT);
                    RD_LOGGER.info("junction canonicalization enabled (max shift: {})",
                            JunctionCanonicalizer.DEFAULT_MAX_SHIFT);
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

    // fails fast on BAM/sidecar mismatch — missing sidecar entries caused silent _tx name leakage into mate fields.
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
                final LiftedMateInfo info = LiftBackRecordOps.toLiftedMateInfo(
                        record, result, mConfig.UnmapAboveNh, mConfig.UnmapBelowMapq);
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

        final LiftBackGroupProcessor processor = new LiftBackGroupProcessor(
                resolver, mJunctionRescueResolver, mSoftclipTailExtender, mTerminalJunctionCollapser, mJunctionCanonicalizer,
                mRefSource, mConfig.UnmapAboveNh, mConfig.UnmapBelowMapq, stats);

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
                final LiftBackGroupProcessor.EmitSink sink = (record, result) ->
                {
                    samWriter.addAlignment(record);
                    try
                    {
                        writer.write(record, result);
                    }
                    catch(IOException e)
                    {
                        throw new RuntimeException(e);
                    }
                };

                final SAMRecordIterator iter = samReader.iterator();
                streamGrouped(iter, processor, liftedMateInfoCache, sink);
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

    // Groups by read name so primary + supps are processed together; rescue needs all split-read components at once.
    private void streamGrouped(
            final SAMRecordIterator iter, final LiftBackGroupProcessor processor,
            final LiftedMateInfoCache liftedMateInfoCache, final LiftBackGroupProcessor.EmitSink sink)
    {
        final List<SAMRecord> buffer = new ArrayList<>();
        String currentName = null;
        while(iter.hasNext())
        {
            final SAMRecord record = iter.next();
            final String name = record.getReadName();
            if(currentName != null && !name.equals(currentName))
            {
                processor.processNameGroup(buffer, liftedMateInfoCache, sink);
                buffer.clear();
            }
            buffer.add(record);
            currentName = name;
        }
        if(!buffer.isEmpty())
        {
            processor.processNameGroup(buffer, liftedMateInfoCache, sink);
        }
    }

    // Synchronized because IndexedFastaSequenceFile is not thread-safe and rescue may run from concurrent partition workers.
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

    public static void main(final String[] args)
    {
        final ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SpliceLiftBackConfig.addConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        new SpliceLiftBack(configBuilder).run();
    }
}
