package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.splice.SpliceCommon.ALT_CONTIG_SUFFIX;

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
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.splice.rescue.AnnotatedJunctionIndex;
import com.hartwig.hmftools.redux.splice.rescue.AnnotatedJunctionLoader;
import com.hartwig.hmftools.redux.splice.rescue.ChrIntron;
import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;
import com.hartwig.hmftools.redux.splice.rescue.RescueRejectReason;
import com.hartwig.hmftools.redux.splice.rescue.RescueStatistics;
import com.hartwig.hmftools.redux.splice.tailextend.TailExtensionStatistics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

// Stage 0 of REDUX when -splice_liftback is set: lifts a bwa-mem2 BAM aligned against ref + transcript
// contigs back to genomic coordinates, producing a coord-sorted+indexed genomic BAM that dedup then
// consumes unchanged. Runs on its OWN producer + worker pool BEFORE createPartitionThreads, the same
// pre-dedup slot as RegionUnmapper -- dedup's genomic-region partition threads are untouched.
//
// This is the concurrent re-shaping of the standalone two-pass SpliceLiftBack.main():
//   - the whole-sample pass-1 LiftedMateInfoCache is DELETED (it was the exp9 OOM). A complete
//     name-sorted BAM keeps both mates + all supps in one contiguous name-group, so a PER-GROUP cache
//     satisfies every partner lookup and the two passes collapse to one. See [[project_redux_migration_scoping]].
//   - parallelism axis is read-NAME chunks (one producer streams contiguous name-groups into a bounded
//     queue, N workers transform), orthogonal to dedup's genomic-region partitions.
//
// Wired into ReduxApplication.run() before setReadLength; no-op unless -splice_liftback is set. The
// orchestration + resource/sort/header/rebind helpers are ported from SpliceLiftBack; what remains is
// the per-group transform (LiftBackWorker.processNameGroup) and the per-worker stat merge.
public class SpliceLiftBackStage
{
    private final ReduxConfig mConfig;
    private final SpliceStageConfig mSpliceConfig;

    // contiguous name-groups are batched into chunks of ~this many reads, cut only at a name boundary.
    static final int CHUNK_TARGET_READS = 5000;

    // bounded in-flight chunk queue so memory scales with concurrency, not sample size (the OOM guard).
    // Capacity is a small multiple of worker count, mirroring MaxCachedReadsPerThread's intent.
    private static final int CHUNK_QUEUE_DEPTH_PER_THREAD = 2;

    // reference-compared sentinel a producer enqueues (one per worker) to signal end-of-stream.
    public static final List<SAMRecord> END_OF_STREAM = new ArrayList<>();

    public SpliceLiftBackStage(final ReduxConfig config, final SpliceStageConfig spliceConfig)
    {
        mConfig = config;
        mSpliceConfig = spliceConfig;
    }

    // entry point for ReduxApplication.run(): no-op + true when the stage is disabled, so REDUX is
    // unchanged. On success rebinds config.BamFiles to the lifted genomic BAM and returns true.
    public static boolean runStage(final ReduxConfig config, final SpliceStageConfig spliceConfig)
    {
        if(spliceConfig == null || !spliceConfig.enabled())
            return true;

        return new SpliceLiftBackStage(config, spliceConfig).run();
    }

    public boolean run()
    {
        if(!mSpliceConfig.enabled())
            return true;

        final long startTimeMs = System.currentTimeMillis();
        RD_LOGGER.info("splice liftback stage starting");

        // 1. name-sort the bwa BAM so primary + mate + supps are contiguous in the stream.
        final String nameSortedBam = nameSortInput();

        // 2. load shared read-only resources once (sidecar, exon/junction indices, ref source,
        //    rescue/extend/collapse engines). Annotation-sized, sample-independent, shared across workers.
        final LiftBackResources resources = buildResources(nameSortedBam);

        // output header derived from the name-sorted BAM with the _tx alt contigs stripped from @SQ.
        final SAMFileHeader outputHeader = buildOutputHeader(nameSortedBam);

        final int workerCount = Math.max(mConfig.Threads, 1);
        final BlockingQueue<List<SAMRecord>> chunkQueue =
                new ArrayBlockingQueue<>(Math.max(workerCount * CHUNK_QUEUE_DEPTH_PER_THREAD, 2));

        // one producer, N workers. Each worker writes its own shard and owns its engines + ref-source handle.
        final ChunkProducer producer = new ChunkProducer(
                nameSortedBam, mConfig.RefGenomeFile, chunkQueue, workerCount);

        final List<LiftBackWorker> workers = Lists.newArrayList();
        final List<Thread> threadTasks = Lists.newArrayList();
        final List<String> shardBams = Lists.newArrayList();
        threadTasks.add(producer);

        for(int i = 0; i < workerCount; ++i)
        {
            final String shardBam = formShardBamPath(i);
            shardBams.add(shardBam);
            final LiftBackWorker worker = new LiftBackWorker(chunkQueue, resources, outputHeader, shardBam);
            workers.add(worker);
            threadTasks.add(worker);
        }

        if(!runThreadTasks(threadTasks))
            return false;

        mergeWorkerStats(workers);

        final String unsortedBam = formUnsortedBamPath();
        if(!concatenateShards(shardBams, unsortedBam))
            return false;

        final String liftedBam = formLiftedBamPath();
        if(!sortAndIndex(unsortedBam, liftedBam))
            return false;

        rebindBamFiles(liftedBam);

        cleanupIntermediates(nameSortedBam, unsortedBam, shardBams);

        RD_LOGGER.info("splice liftback stage complete, mins({})", runTimeMinsStr(startTimeMs));
        return true;
    }

    // shared read-only across workers: the resolver and junction index hold no mutable state. Each worker
    // builds its own engines + ref-source handle from this (RefSequenceSource and the rescue/tail-extend
    // stats are not thread-safe).
    public static final class LiftBackResources
    {
        public final LiftBackResolver Resolver;
        public final AnnotatedJunctionIndex JunctionIndex; // nullable
        public final String RefGenomeFile;
        public final boolean RescueViaSupp;
        public final boolean ExtendSoftclipTails;
        public final int TerminalAnchor;
        public final int UnmapAboveNh;
        public final int UnmapBelowMapq;

        public LiftBackResources(
                final LiftBackResolver resolver, final AnnotatedJunctionIndex junctionIndex, final String refGenomeFile,
                final boolean rescueViaSupp, final boolean extendSoftclipTails, final int terminalAnchor,
                final int unmapAboveNh, final int unmapBelowMapq)
        {
            Resolver = resolver;
            JunctionIndex = junctionIndex;
            RefGenomeFile = refGenomeFile;
            RescueViaSupp = rescueViaSupp;
            ExtendSoftclipTails = extendSoftclipTails;
            TerminalAnchor = terminalAnchor;
            UnmapAboveNh = unmapAboveNh;
            UnmapBelowMapq = unmapBelowMapq;
        }

        public RefSequenceSource openRefSource()
        {
            return SpliceLiftBackStage.openRefSource(RefGenomeFile);
        }
    }

    // name-sort the bwa BAM so primary + mate + supps are contiguous in the stream. Splice liftback
    // takes a single bwa input BAM (not REDUX's multi-BAM merge case).
    private String nameSortInput()
    {
        if(mConfig.BamToolPath == null)
            throw new IllegalStateException("splice liftback name-sort requires -" + BamToolName.BAMTOOL_PATH);

        if(mConfig.BamFiles.size() != 1)
            throw new IllegalStateException("splice liftback expects a single input BAM, got " + mConfig.BamFiles.size());

        final String nameSortedBam = mConfig.OutputDir + "splice_liftback.name_sorted_input.bam";
        RD_LOGGER.info("name-sorting input via {}: {}", fromPath(mConfig.BamToolPath), nameSortedBam);

        final long startTimeMs = System.currentTimeMillis();
        final List<String> command = new ArrayList<>();
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
        command.add(mConfig.BamFiles.get(0));
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


    private LiftBackResources buildResources(final String nameSortedBam)
    {
        final List<ContigEntry> contigEntries = ContigSidecar.read(mSpliceConfig.ContigSidecarFile);
        RD_LOGGER.info("loaded {} contig entries", contigEntries.size());

        ExonRegionIndex exonIndex = null;
        if(mSpliceConfig.EnsemblDataDir != null)
        {
            try
            {
                exonIndex = ExonRegionIndex.load(mSpliceConfig.EnsemblDataDir);
                RD_LOGGER.info("loaded annotated-exon index from {}", mSpliceConfig.EnsemblDataDir);
            }
            catch(IOException e)
            {
                throw new IllegalStateException("failed to load annotated exons: " + e, e);
            }
        }

        final LiftBackResolver resolver = new LiftBackResolver(contigEntries, exonIndex);
        validateBamAgainstSidecar(nameSortedBam, resolver.contigNames());

        AnnotatedJunctionIndex junctionIndex = null;
        if(mSpliceConfig.RescueViaSupp || mSpliceConfig.ExtendSoftclipTails)
        {
            try
            {
                final Set<ChrIntron> annotated = AnnotatedJunctionLoader.load(mSpliceConfig.EnsemblDataDir);
                junctionIndex = new AnnotatedJunctionIndex(annotated);
                RD_LOGGER.info("loaded {} annotated junctions from {}", annotated.size(), mSpliceConfig.EnsemblDataDir);
            }
            catch(IOException e)
            {
                throw new IllegalStateException("failed to load annotated junctions: " + e, e);
            }
        }

        return new LiftBackResources(
                resolver, junctionIndex, mConfig.RefGenomeFile,
                mSpliceConfig.RescueViaSupp, mSpliceConfig.ExtendSoftclipTails, SpliceCommon.MIN_JUNCTION_ANCHOR,
                mSpliceConfig.UnmapAboveNh, mSpliceConfig.UnmapBelowMapq);
    }

    // fails fast on a BAM/sidecar mismatch -- root cause of a prior silent failure where alt contigs
    // missing from the sidecar fell through as ref pass-through and leaked _tx names into mate fields.
    private void validateBamAgainstSidecar(final String nameSortedBam, final Set<String> sidecarContigs)
    {
        final Set<String> missing = new TreeSet<>();
        try(SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(nameSortedBam)))
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
            throw new RuntimeException("failed to read BAM header: " + nameSortedBam, e);
        }

        if(!missing.isEmpty())
        {
            throw new IllegalStateException(String.format(
                    "BAM/sidecar mismatch: %d alt contig(s) in BAM @SQ are absent from sidecar %s -- first few: %s",
                    missing.size(), mSpliceConfig.ContigSidecarFile,
                    missing.stream().limit(5).collect(Collectors.joining(","))));
        }
    }

    // one handle per caller. IndexedFastaSequenceFile is not thread-safe, so workers each open their own.
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

    // clone the name-sorted BAM header, mark it unsorted, and strip the _tx alt contigs from @SQ so the
    // lifted BAM carries a pure genomic dictionary (SpliceLiftBack.java:323-326 + stripAltContigsFromHeader
    // 918-934). RISK: confirm the resulting SQ dictionary matches the ref dict REDUX expects (memory risk #1).
    private SAMFileHeader buildOutputHeader(final String nameSortedBam)
    {
        try(SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(nameSortedBam)))
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
            RD_LOGGER.info("dropped {} alt contig @SQ entries from output header", dropped);
            return header;
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read header from " + nameSortedBam, e);
        }
    }

    private boolean sortAndIndex(final String unsortedBam, final String sortedBam)
    {
        if(mConfig.BamToolPath == null)
        {
            RD_LOGGER.info("no -{} configured; leaving unsorted BAM at {}", BamToolName.BAMTOOL_PATH, unsortedBam);
            return false;
        }

        final BamToolName toolName = fromPath(mConfig.BamToolPath);
        RD_LOGGER.info("sorting BAM via {}: {}", toolName, sortedBam);

        if(!BamOperations.sortBam(toolName, mConfig.BamToolPath, unsortedBam, sortedBam, mConfig.Threads))
            return false;

        if(toolName == BamToolName.SAMTOOLS)
        {
            if(!BamOperations.indexBam(toolName, mConfig.BamToolPath, sortedBam, mConfig.Threads))
                return false;
        }

        return true;
    }

    // mutate config.BamFiles in place so every downstream stage (dedup, unmap, partition readers) reads
    // the lifted genomic BAM. BamFiles is a mutable ArrayList (Collectors.toList in ReduxConfig ctor).
    private void rebindBamFiles(final String liftedBam)
    {
        mConfig.BamFiles.clear();
        mConfig.BamFiles.add(liftedBam);
        RD_LOGGER.info("rebound input BAM to lifted genomic BAM: {}", liftedBam);
    }

    private void mergeWorkerStats(final List<LiftBackWorker> workers)
    {
        final LiftBackStats combined = new LiftBackStats();
        for(final LiftBackWorker worker : workers)
            combined.merge(worker.liftBackStats());
        combined.logSummary();

        if(mSpliceConfig.RescueViaSupp)
            logRescueStats(workers);
        if(mSpliceConfig.ExtendSoftclipTails)
            logTailExtensionStats(workers);

        long collapsedLeading = 0;
        long collapsedTrailing = 0;
        for(final LiftBackWorker worker : workers)
        {
            collapsedLeading += worker.collapsedLeading();
            collapsedTrailing += worker.collapsedTrailing();
        }
        if(collapsedLeading > 0 || collapsedTrailing > 0)
            RD_LOGGER.info("terminal micro-junction collapse: leading={} trailing={}", collapsedLeading, collapsedTrailing);
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
        RD_LOGGER.info("rescue-via-supp summary: candidates={} merged={} suppClamped={} (depth1={} depth2={} depth3={} depth4={})",
                candidates, merged, clamp, depth[1], depth[2], depth[3], depth[4]);
        for(final RescueRejectReason reason : RescueRejectReason.values())
        {
            if(rejects.getOrDefault(reason, 0) > 0)
                RD_LOGGER.info("rescue-via-supp reject {}: {}", reason, rejects.get(reason));
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
        RD_LOGGER.info("extend-softclip-tails summary: evaluated={} extended={} basesLead={} basesTrail={}",
                evaluated, extended, basesLead, basesTrail);
    }

    private boolean concatenateShards(final List<String> shardBams, final String unsortedBam)
    {
        final BamToolName toolName = fromPath(mConfig.BamToolPath);
        RD_LOGGER.info("concatenating {} liftback shards via {}", shardBams.size(), toolName);
        // samtools cat is a plain BGZF block concat and rejects -@, so pass threads=1 (matches FileWriterCache/FinalBamWriter)
        return BamOperations.concatenateBams(toolName, mConfig.BamToolPath, unsortedBam, shardBams, 1);
    }

    private String formShardBamPath(final int index)
    {
        return mConfig.OutputDir + "splice_liftback.shard_" + index + ".bam";
    }

    private String formUnsortedBamPath()
    {
        return mConfig.OutputDir + "splice_liftback.unsorted.bam";
    }

    private String formLiftedBamPath()
    {
        return mConfig.OutputDir + "splice_liftback.lifted.bam";
    }

    private void cleanupIntermediates(final String nameSortedBam, final String unsortedBam, final List<String> shardBams)
    {
        if(mConfig.KeepInterimBams)
            return;

        deleteQuietly(nameSortedBam);
        deleteQuietly(unsortedBam);
        shardBams.forEach(SpliceLiftBackStage::deleteQuietly);
    }

    private static void deleteQuietly(final String path)
    {
        try
        {
            Files.deleteIfExists(Paths.get(path));
        }
        catch(IOException e)
        {
            RD_LOGGER.warn("could not delete intermediate {}: {}", path, e.toString());
        }
    }

    // Single producer: streams the name-sorted BAM and batches contiguous read-name groups into chunks
    // of ~CHUNK_TARGET_READS reads, cutting only at a name boundary so every chunk holds whole name-groups
    // (the per-group cache invariant). Pushes chunks onto the bounded queue, then enqueues one
    // END_OF_STREAM sentinel per worker. Mirrors SpliceLiftBack.streamGrouped's grouping, but never holds
    // more than one chunk + the in-flight queue in memory.
    private static final class ChunkProducer extends Thread
    {
        private final String mNameSortedBam;
        private final String mRefGenomeFile;
        private final BlockingQueue<List<SAMRecord>> mQueue;
        private final int mWorkerCount;

        private ChunkProducer(
                final String nameSortedBam, final String refGenomeFile,
                final BlockingQueue<List<SAMRecord>> queue, final int workerCount)
        {
            mNameSortedBam = nameSortedBam;
            mRefGenomeFile = refGenomeFile;
            mQueue = queue;
            mWorkerCount = workerCount;
        }

        @Override
        public void run()
        {
            try(SamReader samReader = SamReaderFactory.makeDefault()
                    .referenceSequence(new File(mRefGenomeFile))
                    .open(new File(mNameSortedBam)))
            {
                List<SAMRecord> chunk = new ArrayList<>();
                String currentName = null;

                final SAMRecordIterator iter = samReader.iterator();
                while(iter.hasNext())
                {
                    final SAMRecord record = iter.next();
                    final String name = record.getReadName();

                    // at a name boundary the chunk holds only complete groups, so it's safe to cut here
                    // once it's reached the target size.
                    if(currentName != null && !name.equals(currentName) && chunk.size() >= CHUNK_TARGET_READS)
                    {
                        mQueue.put(chunk);
                        chunk = new ArrayList<>();
                    }

                    chunk.add(record);
                    currentName = name;
                }

                if(!chunk.isEmpty())
                    mQueue.put(chunk);

                for(int i = 0; i < mWorkerCount; ++i)
                    mQueue.put(END_OF_STREAM);
            }
            catch(IOException | InterruptedException e)
            {
                RD_LOGGER.error("splice liftback chunk producer failed: {}", e.toString());
                System.exit(1);
            }
        }
    }
}
