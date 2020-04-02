package com.hartwig.hmftools.sage.quality;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class QualityApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(QualityApplication.class);

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException {
        try (final QualityApplication app = new QualityApplication()) {
            app.run();

        } catch (Exception e) {
            System.out.println(e);
        }
    }

    private final ExecutorService executorService;
    private final ReferenceSequenceFile referenceSequenceFile;
    private final List<Future<Collection<QualityCount>>> counters = Lists.newArrayList();

    public QualityApplication() throws IOException {
        this.referenceSequenceFile =
                new IndexedFastaSequenceFile(new File("/Users/jon/hmf/resources/Homo_sapiens.GRCh37.GATK.illumina.fasta"));
        this.executorService = Executors.newFixedThreadPool(7);
    }

    public void run() throws ExecutionException, InterruptedException, IOException {
        LOGGER.info("Starting");

        addAllRegions("17", 50_000_000, 60_000_000);

        final List<QualityCount> allCounts = Lists.newArrayList();
        for (Future<Collection<QualityCount>> counter : counters) {
            allCounts.addAll(counter.get());
        }

        final Collection<QualityCount> simple = QualityGrouping.groupWithoutPosition(allCounts);
        final Collection<QualityCount> qualityCounts = QualityGrouping.groupByQuality(allCounts);
        final Collection<QualityCount> strandCounts = QualityGrouping.groupByStrand(allCounts);
        final Collection<QualityCount> strandOnly = QualityGrouping.groupByStrandOnly(allCounts);

        LOGGER.info("Finishing");
        QualityFile.write("/Users/jon/hmf/analysis/sageValidation/quality/detailed.csv", simple);
//        QualityFile.write("/Users/jon/hmf/analysis/sageValidation/quality/bad_detailed.csv", simple);
//        QualityFile.write("/Users/jon/hmf/analysis/sageValidation/quality/bad_quality.csv", qualityCounts);
//        QualityFile.write("/Users/jon/hmf/analysis/sageValidation/quality/bad_strand.csv", strandCounts);
//        QualityFile.write("/Users/jon/hmf/analysis/sageValidation/quality/bad_strandOnly.csv", strandOnly);

        System.out.println("sdf");

    }

    public void addAllRegions(String contig) {
        addAllRegions(contig, 1, referenceSequenceFile.getSequence(contig).length());
    }

    public void addAllRegions(String contig, int minPosition, int maxPosition) {
        final int regionSliceSize = 500_000;
        for (int i = 0; ; i++) {
            int start = 1 + i * regionSliceSize;
            int end = start + regionSliceSize - 1;

            if (end < minPosition) {
                continue;
            }

            addRegion(contig, Math.max(start, minPosition), Math.min(end, maxPosition));

            if (end >= maxPosition) {
                break;
            }
        }
    }

    @NotNull
    public void addRegion(String contig, int start, int end) {
        final GenomeRegion bounds = GenomeRegions.create(contig, start, end);
//        LOGGER.info("Adding region {}", bounds);

        Future<Collection<QualityCount>> future = executorService.submit(() -> new QualityRegion(
                "/Users/jon/hmf/analysis/COLO829T/bams/COLO829T.chr17.bam",
//                "/Users/jon/hmf/analysis/sageValidation/quality/CPCT02010323T.qual.10M.bam",
                referenceSequenceFile).regionCount(bounds));

        counters.add(future);
    }

    @Override
    public void close() throws Exception {
        executorService.shutdown();
        referenceSequenceFile.close();
    }
}
