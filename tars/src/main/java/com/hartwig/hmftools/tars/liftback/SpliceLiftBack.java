package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.tars.liftback.SpliceLiftBackConfig.SORT_BAMTOOL_PATH;
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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
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
import com.hartwig.hmftools.common.region.ChrBaseRegion;
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
// sorted + indexed. The per-record transform lives in LiftBackGroupProcessor (one per worker, no shared state).
public class SpliceLiftBack
{
    private final SpliceLiftBackConfig mConfig;

    // bwa's name-grouped reads are batched into chunks of ~this many reads before handing to a worker.
    static final int CHUNK_TARGET_READS = 5000;

    // bounded in-flight chunk queue so memory scales with concurrency, not sample size (the OOM guard).
    private static final int CHUNK_QUEUE_DEPTH_PER_THREAD = 2;

    // upper bound on parallel input-reader shards; a few parse threads outpace the workers, more just block.
    private static final int READER_SHARD_CAP = 8;

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

        TARS_LOGGER.info("liftback processing complete, mins({}); sorting + merging shards", runTimeMinsStr(startTimeMs));

        // each shard spans the genome in random coordinate order, so sort the (small) shards in parallel -- each
        // fits in memory with little/no temp -- then stream-merge the sorted shards. Avoids a full-size cat copy
        // and a monolithic sort that spills the whole sample to temp.
        final List<String> sortedShards = sortShards(shardBams);
        if(sortedShards == null)
            throw new RuntimeException("failed to sort liftback BAM shards");

        if(!mergeShards(sortedShards, mConfig.formOutputBam()))
            throw new RuntimeException("failed to merge + index lifted BAM");

        if(mConfig.WriteLiftbackTsv)
        {
            concatenateTsvShards(tsvAShards, LiftBackWriter.TSV_A_HEADER_LINE, mConfig.formTsvAFile());
            concatenateTsvShards(tsvBShards, LiftBackWriter.TSV_B_HEADER_LINE, mConfig.formTsvBFile());
        }

        cleanupIntermediates(shardBams, sortedShards, tsvAShards, tsvBShards);

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

        // a single-thread BGZF parse starves the workers, so the input is read by several shard threads that
        // each parse their own byte range (split on read-name-group boundaries). A handful saturates the queue;
        // beyond that the bounded queue just blocks them, so the count is capped well below the worker pool.
        final int shardCount = Math.max(1, Math.min(workerCount, READER_SHARD_CAP));
        final ShardedChunkProducer producer = new ShardedChunkProducer(
                mConfig.InputBam, mConfig.RefGenomeFile, chunkQueue, workerCount, CHUNK_TARGET_READS, shardCount);

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
            exonIndex = ExonRegionIndex.load(mConfig.EnsemblDataDir, mConfig.RefGenVersion);
            TARS_LOGGER.info("loaded annotated-exon index from {}", mConfig.EnsemblDataDir);
        }

        final LiftBackResolver resolver = new LiftBackResolver(contigEntries, exonIndex);
        validateBamAgainstSidecar(resolver.contigNames());

        AnnotatedJunctionIndex junctionIndex = null;
        if(mConfig.RescueViaSupp || mConfig.ExtendSoftclipTails)
        {
            final Set<ChrBaseRegion> annotated = AnnotatedJunctionLoader.load(mConfig.EnsemblDataDir, mConfig.RefGenVersion);
            junctionIndex = new AnnotatedJunctionIndex(annotated);
            TARS_LOGGER.info("loaded {} annotated junctions from {}", annotated.size(), mConfig.EnsemblDataDir);
        }

        ExcludedRegions excludedRegions = null;
        if(mConfig.RnaUnmapRegionsFile != null)
        {
            excludedRegions = ExcludedRegions.load(mConfig.RnaUnmapRegionsFile);
            TARS_LOGGER.info("loaded excluded regions from {}", mConfig.RnaUnmapRegionsFile);
        }

        return new LiftBackResources(
                resolver, junctionIndex, mConfig.RefGenomeFile,
                mConfig.RescueViaSupp, mConfig.ExtendSoftclipTails, SpliceCommon.MIN_JUNCTION_ANCHOR,
                excludedRegions);
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

        long junctionsCanonicalized = 0;
        long overCapUnmapped = 0;
        for(final LiftBackWorker worker : workers)
        {
            junctionsCanonicalized += worker.junctionsCanonicalized();
            overCapUnmapped += worker.overCapUnmapped();
        }
        TARS_LOGGER.info("junction-canonicalize summary: shifted={}", junctionsCanonicalized);
        if(overCapUnmapped > 0)
            TARS_LOGGER.info("over-cap unmap: {} primaries unmapped (MAPQ 0 + no XA, maps past the bwa XA cap)", overCapUnmapped);

        long excludedReads = 0;
        for(final LiftBackWorker worker : workers)
            excludedReads += worker.excludedReads();
        if(excludedReads > 0)
            TARS_LOGGER.info("excluded-region reads (primary unmapped / supp dropped, post-lift): {}", excludedReads);
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
                TARS_LOGGER.debug("rescue-via-supp reject {}: {}", reason, rejects.get(reason));
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

    // sort each shard concurrently as a subprocess. Each shard is ~1/workerCount of the data, so with a modest
    // per-shard -m it sorts mostly in memory; running workerCount of them at once gives sort parallelism the
    // single-process sort lacks. Peak RAM ~= Threads x sort_memory_gb, so size that flag with the box in mind.
    private List<String> sortShards(final List<String> shardBams)
    {
        if(mConfig.SortBamToolPath == null)
        {
            TARS_LOGGER.error("no -{} or -{} configured; cannot sort lifted shards", BamToolName.BAMTOOL_PATH, SORT_BAMTOOL_PATH);
            return null;
        }

        final BamToolName toolName = fromPath(mConfig.SortBamToolPath);
        final Integer sortMemoryGb = mConfig.SortMemoryGb > 0 ? mConfig.SortMemoryGb : null;
        TARS_LOGGER.info("sorting {} shards via {} ({} concurrent)", shardBams.size(), toolName, mConfig.Threads);

        final List<String> sortedShards = Lists.newArrayList();
        for(int i = 0; i < shardBams.size(); ++i)
            sortedShards.add(formSortedShardPath(i));

        final ExecutorService executor = Executors.newFixedThreadPool(mConfig.Threads);
        final List<Future<Boolean>> futures = Lists.newArrayList();
        for(int i = 0; i < shardBams.size(); ++i)
        {
            final String inputShard = shardBams.get(i);
            final String sortedShard = sortedShards.get(i);
            futures.add(executor.submit(() ->
                    BamOperations.sortBam(toolName, mConfig.SortBamToolPath, inputShard, sortedShard, 1, sortMemoryGb)));
        }
        executor.shutdown();

        boolean allSorted = true;
        for(final Future<Boolean> future : futures)
        {
            try
            {
                if(!future.get())
                    allSorted = false;
            }
            catch(InterruptedException | ExecutionException e)
            {
                TARS_LOGGER.error("shard sort failed: {}", e.toString());
                allSorted = false;
            }
        }

        return allSorted ? sortedShards : null;
    }

    private boolean mergeShards(final List<String> sortedShards, final String outputBam)
    {
        // merge over the coordinate-sorted shards is a streaming k-way merge: near-zero temp, no decode-resort.
        final BamToolName toolName = fromPath(mConfig.BamToolPath);
        TARS_LOGGER.info("merging {} sorted shards via {}: {}", sortedShards.size(), toolName, outputBam);

        if(!BamOperations.mergeBams(toolName, mConfig.BamToolPath, outputBam, sortedShards, mConfig.Threads))
            return false;

        return BamOperations.indexBam(toolName, mConfig.BamToolPath, outputBam, mConfig.Threads);
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
        return mConfig.OutputDir + mConfig.prefix() + ".shard_" + index + ".bam";
    }

    private String formSortedShardPath(final int index)
    {
        return mConfig.OutputDir + mConfig.prefix() + ".shard_" + index + ".sorted.bam";
    }

    private String formTsvShardPath(final int index, final String kind)
    {
        return mConfig.OutputDir + mConfig.prefix() + "." + kind + ".shard_" + index + ".tsv";
    }

    private void cleanupIntermediates(
            final List<String> shardBams, final List<String> sortedShards,
            final List<String> tsvAShards, final List<String> tsvBShards)
    {
        shardBams.forEach(SpliceLiftBack::deleteQuietly);
        sortedShards.forEach(SpliceLiftBack::deleteQuietly);
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
