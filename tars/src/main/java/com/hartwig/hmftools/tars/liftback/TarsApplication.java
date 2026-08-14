package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.tars.common.TarsConstants.ALT_CONTIG_SUFFIX;
import static com.hartwig.hmftools.tars.common.TarsConstants.APP_NAME;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.common.ContigSidecar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReaderFactory;

// Driver for the liftback stage. Cuts bwa's name-grouped BAM (mates + supplementaries contiguous) into whole-fragment
// chunks, lifts them across N workers, then sorts + indexes the concatenated output.
public class TarsApplication
{
    private final TarsConfig mConfig;

    private static final int CHUNK_TARGET_READS = 5000;
    private static final int CHUNK_QUEUE_DEPTH_PER_THREAD = 2;
    private static final int READER_SHARD_CAP = 8;
    private static final double MAX_LIFT_FAILURE_FRACTION = 0.01;

    public TarsApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new TarsConfig(configBuilder);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        List<ContigEntry> contigEntries = ContigSidecar.read(mConfig.ContigSidecarFile);
        SAMFileHeader inputHeader = readInputHeader();

        LiftBackResources resources = buildResources(contigEntries, inputHeader);
        SAMFileHeader outputHeader = buildOutputHeader(inputHeader);

        List<String> shardBams = new ArrayList<>();
        List<LiftBackWorker> workers = runChunkStream(resources, outputHeader, shardBams);
        if(workers == null)
        {
            TARS_LOGGER.error("liftback chunk stream failed");
            System.exit(1);
        }

        RunCounts counts = aggregate(workers);
        logLiftFailureRate(counts);
        writeSummary(counts);

        TARS_LOGGER.info("liftback processing complete, mins({}); concatenating + sorting shards", runTimeMinsStr(startTimeMs));

        String unsortedBam = mConfig.formUnsortedBam();
        if(!concatenateShards(shardBams, unsortedBam))
        {
            TARS_LOGGER.error("failed to concatenate liftback BAM shards");
            System.exit(1);
        }

        if(!sortAndIndex(unsortedBam, mConfig.formOutputBam()))
        {
            TARS_LOGGER.error("failed to sort + index lifted BAM");
            System.exit(1);
        }

        cleanupIntermediates(unsortedBam, shardBams);

        TARS_LOGGER.info("TarsApplication complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    // one producer streams bwa's name-grouped BAM into whole-fragment chunks; N workers lift + emit per-shard
    private List<LiftBackWorker> runChunkStream(
            final LiftBackResources resources, final SAMFileHeader outputHeader, final List<String> shardBams)
    {
        int workerCount = Math.max(mConfig.Threads, 1);
        BlockingQueue<List<SAMRecord>> chunkQueue =
                new ArrayBlockingQueue<>(Math.max(workerCount * CHUNK_QUEUE_DEPTH_PER_THREAD, 2));

        // a single-thread BGZF parse starves the workers, so shard threads each parse their own byte range, split on
        // read-name-group boundaries; a handful saturates the bounded queue, hence the low cap.
        int shardCount = Math.max(1, Math.min(workerCount, READER_SHARD_CAP));
        ShardedChunkProducer producer = new ShardedChunkProducer(
                mConfig.InputBam, mConfig.RefGenomeFile, chunkQueue, workerCount, CHUNK_TARGET_READS, shardCount);

        List<LiftBackWorker> workers = new ArrayList<>();
        List<Thread> threadTasks = new ArrayList<>();
        threadTasks.add(producer);

        for(int i = 0; i < workerCount; ++i)
        {
            String shardBam = formShardBamPath(i);
            shardBams.add(shardBam);

            LiftBackWorker worker = new LiftBackWorker(chunkQueue, resources, outputHeader, shardBam);
            workers.add(worker);
            threadTasks.add(worker);
        }

        return runThreadTasks(threadTasks) ? workers : null;
    }

    private LiftBackResources buildResources(final List<ContigEntry> contigEntries, final SAMFileHeader inputHeader)
    {
        // exon + junction annotation are derived from the sidecar's exonSpans, so liftback needs no ensembl cache.
        EnsemblAnnotationIndex annotationIndex = EnsemblAnnotationIndex.fromContigEntries(contigEntries);
        TARS_LOGGER.info("built annotation index from sidecar: {} junctions", annotationIndex.junctionCount());

        // annotation-only rows have no contig to lift against, so the discriminator sees only real contig entries.
        List<ContigEntry> liftEntries = contigEntries.stream()
                .filter(entry -> entry.contigStart() > 0).collect(Collectors.toList());
        LiftBackDiscriminator discriminator = new LiftBackDiscriminator(liftEntries, annotationIndex);
        validateBamAgainstSidecar(inputHeader, discriminator.contigTranslator().contigNames());

        ExcludedRegions excludedRegions = null;
        if(mConfig.RnaUnmapRegionsFile != null)
        {
            excludedRegions = ExcludedRegions.load(mConfig.RnaUnmapRegionsFile);
            TARS_LOGGER.info("loaded excluded regions from {}", mConfig.RnaUnmapRegionsFile);
        }

        return new LiftBackResources(
                discriminator, annotationIndex, mConfig.RefGenomeFile,
                mConfig.Supplementary, excludedRegions);
    }

    // an alt contig missing from the sidecar cannot be lifted and would leak its _tx name into the output, so fail fast
    private void validateBamAgainstSidecar(final SAMFileHeader inputHeader, final Set<String> sidecarContigs)
    {
        Set<String> missing = new TreeSet<>();
        for(SAMSequenceRecord sequenceRecord : inputHeader.getSequenceDictionary().getSequences())
        {
            String name = sequenceRecord.getSequenceName();
            if(name.endsWith(ALT_CONTIG_SUFFIX) && !sidecarContigs.contains(name))
            {
                missing.add(name);
            }
        }

        if(missing.isEmpty())
        {
            return;
        }

        TARS_LOGGER.error(
                "BAM/sidecar mismatch: {} alt contig(s) in BAM @SQ are absent from sidecar {} - first few: {}",
                missing.size(), mConfig.ContigSidecarFile,
                missing.stream().limit(5).collect(Collectors.joining(",")));
        System.exit(1);
    }

    // input @SQ, read once: both the sidecar check and the output header derive from it
    private SAMFileHeader readInputHeader()
    {
        return SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.InputBam))
                .getFileHeader();
    }

    // strip the _tx alt contigs from @SQ so the lifted BAM carries a pure genomic dictionary
    private static SAMFileHeader buildOutputHeader(final SAMFileHeader inputHeader)
    {
        SAMFileHeader header = inputHeader.clone();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        List<SAMSequenceRecord> kept = new ArrayList<>();
        int dropped = 0;
        for(SAMSequenceRecord sequenceRecord : header.getSequenceDictionary().getSequences())
        {
            if(sequenceRecord.getSequenceName().endsWith(ALT_CONTIG_SUFFIX))
            {
                ++dropped;
            }
            else
            {
                kept.add(sequenceRecord);
            }
        }

        header.setSequenceDictionary(new SAMSequenceDictionary(kept));
        TARS_LOGGER.info("dropped {} alt contig @SQ entries from output header", dropped);
        return header;
    }

    // the three deliberate unmapping routes are normal outcomes, kept apart from LiftFailed which alone signals a sidecar/FASTA mismatch
    private record RunCounts(
            long RecordsSeen, long PrimariesSeen, long LiftFailed, long ExcludedRegion, long OverCap, long LowAlignmentScore,
            long SupplementaryCandidates, long PrimaryRevisions, long SupplementaryMerges, long SupplementariesAbsorbed)
    {
    }

    private static RunCounts aggregate(final List<LiftBackWorker> workers)
    {
        long recordsSeen = 0, primariesSeen = 0, liftFailed = 0, excludedRegion = 0, overCap = 0, lowAlignmentScore = 0;
        long suppCandidates = 0, primaryRevisions = 0, suppMerges = 0, suppAbsorbed = 0;

        for(LiftBackWorker worker : workers)
        {
            recordsSeen += worker.recordsSeen();
            primariesSeen += worker.primariesSeen();
            liftFailed += worker.primariesLiftFailed();
            excludedRegion += worker.primariesUnmappedExcludedRegion();
            overCap += worker.primariesUnmappedOverCap();
            lowAlignmentScore += worker.primariesUnmappedLowAlignmentScore();
            suppCandidates += worker.supplementaryCandidates();
            primaryRevisions += worker.primaryRevisions();
            suppMerges += worker.supplementaryMerges();
            suppAbsorbed += worker.supplementariesAbsorbed();
        }

        return new RunCounts(
                recordsSeen, primariesSeen, liftFailed, excludedRegion, overCap, lowAlignmentScore,
                suppCandidates, primaryRevisions, suppMerges, suppAbsorbed);
    }

    // A wholesale lift failure is systemic: almost always a sidecar built against a different FASTA than the reads were
    // aligned to. Deliberate unmapping is excluded, so the rate covers only lifts that produced nothing. Logged, not
    // fatal, so the written BAM stays inspectable.
    private static void logLiftFailureRate(final RunCounts counts)
    {
        TARS_LOGGER.info(
                "processed {} records, {} mapped primaries, {} failed to lift",
                counts.RecordsSeen, counts.PrimariesSeen, counts.LiftFailed);

        if(counts.PrimariesSeen > 0 && counts.LiftFailed / (double) counts.PrimariesSeen > MAX_LIFT_FAILURE_FRACTION)
        {
            TARS_LOGGER.error(
                    "lift failure rate {} / {} = {}% exceeds {}% - likely sidecar/FASTA mismatch",
                    counts.LiftFailed, counts.PrimariesSeen,
                    String.format("%.2f", 100.0 * counts.LiftFailed / counts.PrimariesSeen),
                    String.format("%.2f", 100.0 * MAX_LIFT_FAILURE_FRACTION));
        }
    }

    private void writeSummary(final RunCounts counts)
    {
        String filename = mConfig.formSummaryFile();

        try(BufferedWriter writer = createBufferedWriter(filename))
        {
            writer.write(String.join(TSV_DELIM, "Metric", "Value", "Pct", "Basis"));
            writer.newLine();

            writeTotal(writer, "records_total", counts.RecordsSeen);
            writeTotal(writer, "mapped_primaries", counts.PrimariesSeen);

            writeMetric(writer, "lift_failed", counts.LiftFailed, "mapped_primaries", counts.PrimariesSeen);
            writeMetric(writer, "unmapped_excluded_region", counts.ExcludedRegion, "mapped_primaries", counts.PrimariesSeen);
            writeMetric(writer, "unmapped_over_cap", counts.OverCap, "mapped_primaries", counts.PrimariesSeen);
            writeMetric(writer, "unmapped_low_alignment_score", counts.LowAlignmentScore, "mapped_primaries", counts.PrimariesSeen);

            writeMetric(
                    writer, "supp_merge_candidates", counts.SupplementaryCandidates, "mapped_primaries", counts.PrimariesSeen);
            writeMetric(writer, "primary_revised", counts.PrimaryRevisions, "mapped_primaries", counts.PrimariesSeen);
            writeMetric(
                    writer, "supp_merged", counts.SupplementaryMerges, "supp_merge_candidates", counts.SupplementaryCandidates);
            writeMetric(
                    writer, "supps_absorbed", counts.SupplementariesAbsorbed, "supp_merge_candidates",
                    counts.SupplementaryCandidates);

            TARS_LOGGER.info("wrote summary to {}", filename);
        }
        catch(IOException e)
        {
            TARS_LOGGER.warn("failed to write summary {}: {}", filename, e.toString());
        }
    }

    private static void writeTotal(final BufferedWriter writer, final String metric, final long value) throws IOException
    {
        writer.write(String.join(TSV_DELIM, metric, String.valueOf(value), "", ""));
        writer.newLine();
    }

    // Basis names the denominator: unmap reasons are shares of mapped primaries, merge outcomes shares of merge candidates.
    private static void writeMetric(
            final BufferedWriter writer, final String metric, final long value, final String basisName, final long basis)
            throws IOException
    {
        String pct = basis > 0 ? String.format("%.4f", 100.0 * value / basis) : "";
        writer.write(String.join(TSV_DELIM, metric, String.valueOf(value), pct, basisName));
        writer.newLine();
    }

    private boolean concatenateShards(final List<String> shardBams, final String unsortedBam)
    {
        BamToolName toolName = fromPath(mConfig.BamToolPath);
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

        BamToolName toolName = fromPath(mConfig.BamToolPath);
        TARS_LOGGER.info("sorting BAM via {}: {}", toolName, sortedBam);

        if(!BamOperations.sortBam(toolName, mConfig.BamToolPath, unsortedBam, sortedBam, mConfig.Threads))
        {
            return false;
        }

        // sambamba sort indexes inline; only samtools needs the explicit index pass.
        if(toolName == BamToolName.SAMTOOLS)
        {
            if(!BamOperations.indexBam(toolName, mConfig.BamToolPath, sortedBam, mConfig.Threads))
            {
                return false;
            }
        }

        return true;
    }

    private String formShardBamPath(final int index)
    {
        return mConfig.filePrefix() + ".shard_" + index + ".bam";
    }

    private void cleanupIntermediates(final String unsortedBam, final List<String> shardBams)
    {
        deleteQuietly(unsortedBam);
        shardBams.forEach(TarsApplication::deleteQuietly);
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
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        TarsConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        TarsApplication application = new TarsApplication(configBuilder);
        application.run();
    }
}
