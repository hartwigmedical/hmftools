package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.cli.Configs.defaultIntValue;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationRegions;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SageQualityRecalibrationBedApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SageQualityRecalibrationBedApplication.class);

    private static final String BQR_SAMPLE_SIZE = "bqr_sample_size";
    private static final String REF_GENOME = "ref_genome";
    private static final String OUT = "out";
    private static final int DEFAULT_BQR_SAMPLE_SIZE = 2_000_000;

    public static void main(final String[] args) {
        final Options options = createOptions();
        try (final SageQualityRecalibrationBedApplication application = new SageQualityRecalibrationBedApplication(options, args)) {
            application.run();
        } catch (Exception e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageQualityRecalibrationBedApplication", options);
            System.exit(1);
        }
    }

    private final IndexedFastaSequenceFile refGenome;
    private final int sampleSize;
    private final String out;

    private SageQualityRecalibrationBedApplication(final Options options, final String... args)
            throws ParseException, FileNotFoundException {
        final CommandLine cmd = createCommandLine(args, options);

        this.refGenome = new IndexedFastaSequenceFile(new File(cmd.getOptionValue(REF_GENOME)));
        this.sampleSize = defaultIntValue(cmd, BQR_SAMPLE_SIZE, DEFAULT_BQR_SAMPLE_SIZE);
        this.out = cmd.getOptionValue(OUT);
    }

    public void run() throws IOException {
        List<GenomeRegion> regions = new QualityRecalibrationRegions(refGenome).regions(sampleSize);
        NamedBedFile.writeUnnamedBedFile(out, regions);
    }

    @Override
    public void close() throws Exception {
        refGenome.close();
    }

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(REF_GENOME, true, "Path to indexed ref genome fasta file");
        options.addOption(BQR_SAMPLE_SIZE, true, "BQR sampling size per autosome [" + DEFAULT_BQR_SAMPLE_SIZE + "]");

        options.addOption(OUT, true, "Path to output bed");

        return options;
    }

    @NotNull
    static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
