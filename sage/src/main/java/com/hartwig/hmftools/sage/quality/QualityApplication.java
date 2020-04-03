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
import com.hartwig.hmftools.common.cli.Configs;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class QualityApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(QualityApplication.class);

    private static final String OUT = "out";
    private static final String IN = "in";
    private static final String REF_GENOME = "ref_genome";
    private static final String THREADS = "threads";

    public static void main(final String... args) throws ParseException {

        final Options options = new Options();
        options.addOption(OUT, true, "output");
        options.addOption(IN, true, "input");
        options.addOption(REF_GENOME, true, "reference");
        options.addOption(THREADS, true, "reference");

        final CommandLine commandLine = createCommandLine(args, options);

        try (final QualityApplication app = new QualityApplication(commandLine.getOptionValue(REF_GENOME),
                Configs.defaultIntValue(commandLine, THREADS, 32))) {
            app.run(commandLine.getOptionValue(IN), commandLine.getOptionValue(OUT));

        } catch (Exception e) {
            System.out.println(e);
        }
    }

    private final ExecutorService executorService;
    private final ReferenceSequenceFile referenceSequenceFile;
    private final List<Future<Collection<QualityCount>>> counters = Lists.newArrayList();

    public QualityApplication(String referenence, int threads) throws IOException {
        this.referenceSequenceFile = new IndexedFastaSequenceFile(new File(referenence));
        this.executorService = Executors.newFixedThreadPool(threads);
    }

    public void run(String input, String output) throws ExecutionException, InterruptedException, IOException {
        long time = System.currentTimeMillis();

        LOGGER.info("Starting");

        for (final SAMSequenceRecord sequenceRecord : referenceSequenceFile.getSequenceDictionary().getSequences()) {
            if (HumanChromosome.contains(sequenceRecord.getSequenceName())) {
                final HumanChromosome chromosome = HumanChromosome.fromString(sequenceRecord.getSequenceName());
                if (chromosome.isAutosome()) {
                    int start = sequenceRecord.getSequenceLength() - 2_000_000;
                    int end = sequenceRecord.getSequenceLength() - 1_000_001;
                    addAllRegions(input, sequenceRecord.getSequenceName(), start, end);
                }
            }
        }

        List<QualityCount> allCounts = Lists.newArrayList();
        for (Future<Collection<QualityCount>> counter : counters) {
            allCounts.addAll(counter.get());
            allCounts = QualityGrouping.groupWithoutPosition(allCounts);
        }

        QualityFile.write(output, allCounts);
        LOGGER.info("Finished in {} seconds", (System.currentTimeMillis() - time) / 1000);
    }

    @NotNull
    public void addRegion(String bam, String contig, int start, int end) {
        final GenomeRegion bounds = GenomeRegions.create(contig, start, end);
        final Future<Collection<QualityCount>> future =
                executorService.submit(() -> new QualityRegion(bam, referenceSequenceFile).regionCount(bounds));
        counters.add(future);
    }

    public void addAllRegions(String bam, String contig, int minPosition, int maxPosition) {
        final int regionSliceSize = 100_000;
        for (int i = 0; ; i++) {
            int start = minPosition + i * regionSliceSize;
            int end = start + regionSliceSize - 1;

            if (end < minPosition) {
                continue;
            }

            addRegion(bam, contig, Math.max(start, minPosition), Math.min(end, maxPosition));

            if (end >= maxPosition) {
                break;
            }
        }
    }

    @Override
    public void close() throws Exception {
        executorService.shutdown();
        referenceSequenceFile.close();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
