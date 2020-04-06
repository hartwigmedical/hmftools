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

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

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
    private final IndexedFastaSequenceFile refGenome;
    private final List<Future<Collection<QualityCounter>>> counters = Lists.newArrayList();

    public QualityApplication(String referenence, int threads) throws IOException {
        this.refGenome = new IndexedFastaSequenceFile(new File(referenence));
        this.executorService = Executors.newFixedThreadPool(threads);
    }

    public void run(String input, String output) throws ExecutionException, InterruptedException, IOException {
        long time = System.currentTimeMillis();

        LOGGER.info("Starting");

      final QualityRecalibration recalibration = new QualityRecalibration(executorService, refGenome);
      final List<QualityRecalibrationRecord> recalibrationRecords = recalibration.qualityRecalibrationRecords(input).join();
        QualityRecalibrationFile.write(output, recalibrationRecords);
        LOGGER.info("Finished in {} seconds", (System.currentTimeMillis() - time) / 1000);
    }

    @Override
    public void close() throws Exception {
        executorService.shutdown();
        refGenome.close();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
