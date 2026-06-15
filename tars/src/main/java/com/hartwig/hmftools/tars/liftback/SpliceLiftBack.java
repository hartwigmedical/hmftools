package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.tars.common.SpliceCommon.ALT_CONTIG_SUFFIX;
import static com.hartwig.hmftools.tars.common.TarsConfig.APP_NAME;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.common.ContigSidecar;
import com.hartwig.hmftools.tars.common.SpliceCommon;
import com.hartwig.hmftools.tars.liftback.rescue.AnnotatedJunctionIndex;
import com.hartwig.hmftools.tars.liftback.rescue.AnnotatedJunctionLoader;
import com.hartwig.hmftools.tars.liftback.rescue.ChrIntron;
import com.hartwig.hmftools.tars.liftback.rescue.RescueRejectReason;
import com.hartwig.hmftools.tars.liftback.rescue.RescueStatistics;
import com.hartwig.hmftools.tars.liftback.tailextend.TailExtensionStatistics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// Standalone driver over the bwa-mem2 BAM aligned against ref + tx contigs. Consumes bwa's name-grouped
// output directly (mates + supplementaries already contiguous, FASTQ order): ChunkProducer streams it once
// and cuts whole-fragment chunks at name boundaries, N LiftBackWorkers lift + emit each chunk to its own
// BAM + (headerless) TSV shard. No input sort, no index, no fragment cache; only the lifted OUTPUT is
// sorted + indexed. The per-record transform lives in LiftBackGroupProcessor (shared by every worker).
public class SpliceLiftBack
{
    private final SpliceLiftBackConfig mConfig;

    // bwa's name-grouped reads are batched into chunks of ~this many reads before handing to a worker.
    static final int CHUNK_TARGET_READS = 5000;

    // bounded in-flight chunk queue so memory scales with concurrency, not sample size (the OOM guard).
    private static final int CHUNK_QUEUE_DEPTH_PER_THREAD = 2;

    public SpliceLiftBack(final ConfigBuilder configBuilder)
    {
        mConfig = new SpliceLiftBackConfig(configBuilder);
    }

    public void run()
    {
        final long startTimeMs = System.currentTimeMillis();

        final List<ContigEntry> contigEntries = ContigSidecar.read(mConfig.ContigSidecarFile);
        TARS_LOGGER.info("loaded {} contig entries", contigEntries.size());

        final LiftBackResources resources = buildResources(contigEntries);
        final SAMFileHeader outputHeader = buildOutputHeader();

        final List<String> shardBams = Lists.newArrayList();
        final List<String> tsvAShards = Lists.newArrayList();
        final List<String> tsvBShards = Lists.newArrayList();

        final List<LiftBackWorker> workers = runChunkStream(resources, outputHeader, shardBams, tsvAShards, tsvBShards);
        if(workers == null)
            throw new RuntimeException("liftback chunk stream failed");

        mergeWorkerStats(workers);

        final String unsortedBam = mConfig.formUnsortedBam();
        if(!concatenateShards(shardBams, unsortedBam))
            throw new RuntimeException("failed to concatenate liftback BAM shards");

        if(!sortAndIndex(unsortedBam, mConfig.formOutputBam()))
            throw new RuntimeException("failed to sort + index lifted BAM");

        if(mConfig.WriteLiftbackTsv)
        {
            concatenateTsvShards(tsvAShards, LiftBackWriter.TSV_A_HEADER_LINE, mConfig.formTsvAFile());
            concatenateTsvShards(tsvBShards, LiftBackWriter.TSV_B_HEADER_LINE, mConfig.formTsvBFile());
        }

        cleanupIntermediates(unsortedBam, shardBams, tsvAShards, tsvBShards);

        TARS_LOGGER.info("SpliceLiftBack complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    // One producer streams bwa's name-grouped BAM, cutting whole-fragment chunks; N workers lift + emit each
    // chunk to its own BAM + TSV shard. No input sort, no index, no fragment cache.
    private List<LiftBackWorker> runChunkStream(
            final LiftBackResources resources, final SAMFileHeader outputHeader,
            final List<String> shardBams, final List<String> tsvAShards, final List<String> tsvBShards)
    {
        final int workerCount = Math.max(mConfig.Threads, 1);
        final BlockingQueue<List<SAMRecord>> chunkQueue =
                new ArrayBlockingQueue<>(Math.max(workerCount * CHUNK_QUEUE_DEPTH_PER_THREAD, 2));

        // decompress the input across a few of the bam tool's threads so the single producer doesn't starve
        // the workers; the workers still get the bulk of the cores.
        final int decompressThreads = Math.max(1, Math.min(4, mConfig.Threads));
        final ChunkProducer producer = new ChunkProducer(
                mConfig.InputBam, mConfig.RefGenomeFile, chunkQueue, workerCount, CHUNK_TARGET_READS,
                mConfig.BamToolPath, decompressThreads);

        final List<LiftBackWorker> workers = Lists.newArrayList();
        final List<Thread> threadTasks = Lists.newArrayList();
        threadTasks.add(producer);

        for(int i = 0; i < workerCount; ++i)
        {
            final String shardBam = formShardBamPath(i);
            shardBams.add(shardBam);

            // null shard paths disable the worker's TSV writer; debug TSVs are opt-in (whole-sample is huge).
            final String tsvAShard = mConfig.WriteLiftbackTsv ? formTsvShardPath(i, "records") : null;
            final String tsvBShard = mConfig.WriteLiftbackTsv ? formTsvShardPath(i, "alignments") : null;
            if(mConfig.WriteLiftbackTsv)
            {
                tsvAShards.add(tsvAShard);
                tsvBShards.add(tsvBShard);
            }

            final LiftBackWorker worker = new LiftBackWorker(chunkQueue, resources, outputHeader, shardBam, tsvAShard, tsvBShard);
            workers.add(worker);
            threadTasks.add(worker);
        }

        return runThreadTasks(threadTasks) ? workers : null;
    }

    private LiftBackResources buildResources(final List<ContigEntry> contigEntries)
    {
        ExonRegionIndex exonIndex = null;
        if(mConfig.EnsemblDataDir != null)
        {
            try
            {
                exonIndex = ExonRegionIndex.load(mConfig.EnsemblDataDir);
                TARS_LOGGER.info("loaded annotated-exon index from {}", mConfig.EnsemblDataDir);
            }
            catch(IOException e)
            {
                throw new IllegalStateException("failed to load annotated exons: " + e, e);
            }
        }

        final LiftBackResolver resolver = new LiftBackResolver(contigEntries, exonIndex);
        validateBamAgainstSidecar(resolver.contigNames());

        AnnotatedJunctionIndex junctionIndex = null;
        if(mConfig.RescueViaSupp || mConfig.ExtendSoftclipTails)
        {
            try
            {
                final Set<ChrIntron> annotated = AnnotatedJunctionLoader.load(mConfig.EnsemblDataDir);
                junctionIndex = new AnnotatedJunctionIndex(annotated);
                TARS_LOGGER.info("loaded {} annotated junctions from {}", annotated.size(), mConfig.EnsemblDataDir);
            }
            catch(IOException e)
            {
                throw new IllegalStateException("failed to load annotated junctions: " + e, e);
            }
        }

        ExcludedRegions excludedRegions = null;
        if(mConfig.RnaUnmapRegionsFile != null)
        {
            excludedRegions = ExcludedRegions.load(mConfig.RnaUnmapRegionsFile);
            TARS_LOGGER.info("loaded pre-liftback excluded regions from {}", mConfig.RnaUnmapRegionsFile);
        }

        return new LiftBackResources(
                resolver, junctionIndex, mConfig.RefGenomeFile,
                mConfig.RescueViaSupp, mConfig.ExtendSoftclipTails, SpliceCommon.MIN_JUNCTION_ANCHOR,
                mConfig.UnmapAboveNh, mConfig.UnmapBelowMapq, excludedRegions);
    }

    // fails fast on a BAM/sidecar mismatch -- root cause of a prior silent failure where alt contigs missing
    // from the sidecar fell through as ref pass-through and leaked _tx names into mate fields.
    private void validateBamAgainstSidecar(final Set<String> sidecarContigs)
    {
        final Set<String> missing = new TreeSet<>();
        try(SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.InputBam)))
        {
            for(final SAMSequenceRecord sq : reader.getFileHeader().getSequenceDictionary().getSequences())
            {
                final String name = sq.getSequenceName();
                if(name.endsWith(ALT_CONTIG_SUFFIX) && !sidecarContigs.contains(name))
                    missing.add(name);
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read BAM header: " + mConfig.InputBam, e);
        }

        if(!missing.isEmpty())
        {
            throw new IllegalStateException(String.format(
                    "BAM/sidecar mismatch: %d alt contig(s) in BAM @SQ are absent from sidecar %s -- first few: %s",
                    missing.size(), mConfig.ContigSidecarFile,
                    missing.stream().limit(5).collect(Collectors.joining(","))));
        }
    }

    // clone the input BAM header, mark it unsorted, and strip the _tx alt contigs from @SQ so the lifted BAM
    // carries a pure genomic dictionary.
    private SAMFileHeader buildOutputHeader()
    {
        try(SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.InputBam)))
        {
            final SAMFileHeader header = reader.getFileHeader().clone();
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);

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
            TARS_LOGGER.info("dropped {} alt contig @SQ entries from output header", dropped);
            return header;
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read header from " + mConfig.InputBam, e);
        }
    }

    private void mergeWorkerStats(final List<LiftBackWorker> workers)
    {
        final LiftBackStats combined = new LiftBackStats();
        for(final LiftBackWorker worker : workers)
            combined.merge(worker.liftBackStats());
        combined.logSummary();

        try
        {
            combined.writeSummary(mConfig.formSummaryFile());
        }
        catch(IOException e)
        {
            TARS_LOGGER.warn("failed to write liftback summary: {}", e.toString());
        }

        if(mConfig.RescueViaSupp)
            logRescueStats(workers);
        if(mConfig.ExtendSoftclipTails)
            logTailExtensionStats(workers);

        long collapsedLeading = 0;
        long collapsedTrailing = 0;
        for(final LiftBackWorker worker : workers)
        {
            collapsedLeading += worker.collapsedLeading();
            collapsedTrailing += worker.collapsedTrailing();
        }
        if(collapsedLeading > 0 || collapsedTrailing > 0)
            TARS_LOGGER.info("terminal micro-junction collapse: leading={} trailing={}", collapsedLeading, collapsedTrailing);

        long excludedReads = 0;
        for(final LiftBackWorker worker : workers)
            excludedReads += worker.excludedReads();
        if(excludedReads > 0)
            TARS_LOGGER.info("pre-liftback excluded-region reads dropped: {}", excludedReads);
    }

    private static void logRescueStats(final List<LiftBackWorker> workers)
    {
        int candidates = 0;
        int merged = 0;
        int clamp = 0;
        final int[] depth = new int[5];
        final EnumMap<RescueRejectReason,Integer> rejects = new EnumMap<>(RescueRejectReason.class);
        for(final LiftBackWorker worker : workers)
        {
            final RescueStatistics stats = worker.rescueStatistics();
            if(stats == null)
                continue;
            candidates += stats.candidatesEvaluated();
            merged += stats.mergedTotal();
            clamp += stats.suppClampApplied();
            for(int d = 1; d <= 4; ++d)
                depth[d] += stats.mergedAtChainDepth(d);
            for(final RescueRejectReason reason : RescueRejectReason.values())
                rejects.merge(reason, stats.rejectCount(reason), Integer::sum);
        }
        TARS_LOGGER.info("rescue-via-supp summary: candidates={} merged={} suppClamped={} (depth1={} depth2={} depth3={} depth4={})",
                candidates, merged, clamp, depth[1], depth[2], depth[3], depth[4]);
        for(final RescueRejectReason reason : RescueRejectReason.values())
        {
            if(rejects.getOrDefault(reason, 0) > 0)
                TARS_LOGGER.info("rescue-via-supp reject {}: {}", reason, rejects.get(reason));
        }
    }

    private static void logTailExtensionStats(final List<LiftBackWorker> workers)
    {
        int evaluated = 0;
        int extended = 0;
        int basesLead = 0;
        int basesTrail = 0;
        for(final LiftBackWorker worker : workers)
        {
            final TailExtensionStatistics stats = worker.tailExtStatistics();
            if(stats == null)
                continue;
            evaluated += stats.recordsEvaluated();
            extended += stats.recordsExtended();
            basesLead += stats.basesExtendedLead();
            basesTrail += stats.basesExtendedTrail();
        }
        TARS_LOGGER.info("extend-softclip-tails summary: evaluated={} extended={} basesLead={} basesTrail={}",
                evaluated, extended, basesLead, basesTrail);
    }

    private boolean concatenateShards(final List<String> shardBams, final String unsortedBam)
    {
        final BamToolName toolName = fromPath(mConfig.BamToolPath);
        TARS_LOGGER.info("concatenating {} liftback shards via {}", shardBams.size(), toolName);
        // samtools cat is a plain BGZF block concat and rejects -@, so pass threads=1.
        return BamOperations.concatenateBams(toolName, mConfig.BamToolPath, unsortedBam, shardBams, 1);
    }

    private boolean sortAndIndex(final String unsortedBam, final String sortedBam)
    {
        if(mConfig.BamToolPath == null)
        {
            TARS_LOGGER.info("no -{} configured; leaving unsorted BAM at {}", BamToolName.BAMTOOL_PATH, unsortedBam);
            return false;
        }

        final BamToolName toolName = fromPath(mConfig.BamToolPath);
        TARS_LOGGER.info("sorting BAM via {}: {}", toolName, sortedBam);

        if(!BamOperations.sortBam(toolName, mConfig.BamToolPath, unsortedBam, sortedBam, mConfig.Threads))
            return false;

        // sambamba sort indexes inline; only samtools needs the explicit index pass.
        if(toolName == BamToolName.SAMTOOLS)
        {
            if(!BamOperations.indexBam(toolName, mConfig.BamToolPath, sortedBam, mConfig.Threads))
                return false;
        }

        return true;
    }

    // write the header line once, then append each worker's headerless data shard.
    static void concatenateTsvShards(final List<String> shards, final String headerLine, final String outPath)
    {
        try(BufferedWriter writer = FileWriterUtils.createBufferedWriter(outPath))
        {
            writer.write(headerLine);
            writer.newLine();
            for(final String shard : shards)
            {
                try(BufferedReader reader = Files.newBufferedReader(Paths.get(shard)))
                {
                    String line;
                    while((line = reader.readLine()) != null)
                    {
                        writer.write(line);
                        writer.newLine();
                    }
                }
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to concatenate TSV shards into " + outPath, e);
        }
    }

    private String formShardBamPath(final int index)
    {
        return mConfig.OutputDir + "splice_liftback.shard_" + index + ".bam";
    }

    private String formTsvShardPath(final int index, final String kind)
    {
        return mConfig.OutputDir + "splice_liftback." + kind + ".shard_" + index + ".tsv";
    }

    private void cleanupIntermediates(
            final String unsortedBam, final List<String> shardBams,
            final List<String> tsvAShards, final List<String> tsvBShards)
    {
        deleteQuietly(unsortedBam);
        shardBams.forEach(SpliceLiftBack::deleteQuietly);
        tsvAShards.forEach(SpliceLiftBack::deleteQuietly);
        tsvBShards.forEach(SpliceLiftBack::deleteQuietly);
    }

    private static void deleteQuietly(final String path)
    {
        try
        {
            Files.deleteIfExists(Paths.get(path));
        }
        catch(IOException e)
        {
            TARS_LOGGER.warn("could not delete intermediate {}: {}", path, e.toString());
        }
    }

    public static void main(final String[] args)
    {
        final ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SpliceLiftBackConfig.addConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        new SpliceLiftBack(configBuilder).run();
    }
}
