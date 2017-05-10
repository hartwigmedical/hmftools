package com.hartwig.hmftools.purple;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.slicing.SlicerFactory.sortedSlicer;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.freec.FreecCopyNumber;
import com.hartwig.hmftools.common.freec.FreecCopyNumberFactory;
import com.hartwig.hmftools.common.freec.FreecFileLoader;
import com.hartwig.hmftools.common.freec.FreecRatio;
import com.hartwig.hmftools.common.freec.FreecRatioFactory;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.purple.EnrichedCopyNumber;
import com.hartwig.hmftools.common.purple.EnrichedCopyNumberFactory;
import com.hartwig.hmftools.common.purple.FittedCopyNumber;
import com.hartwig.hmftools.common.purple.FittedCopyNumberFactory;
import com.hartwig.hmftools.common.purple.FittedCopyNumberWriter;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityFactory;
import com.hartwig.hmftools.common.purple.FittedPurityWriter;
import com.hartwig.hmftools.common.purple.PadCopyNumber;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFGermlineFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PurityPloidyEstimateApplication {

    private static final Logger LOGGER = LogManager.getLogger(PurityPloidyEstimateApplication.class);

    // Constants
    private static final double MIN_REF_ALLELE_FREQUENCY = 0.4;
    private static final double MAX_REF_ALLELE_FREQUENCY = 0.65;
    private static final int MIN_COMBINED_DEPTH = 10;
    private static final int MAX_COMBINED_DEPTH = 100;
    private static final int MAX_PLOIDY = 12;
    private static final double MIN_PURITY = 0.1;
    private static final double MAX_PURITY = 1.0;
    private static final double PURITY_INCREMENTS = 0.01;
    private static final double MIN_NORM_FACTOR = 0.33;
    private static final double MAX_NORM_FACTOR = 2;
    private static final double NORM_FACTOR_INCREMENTS = 0.01;

    // Options
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String BED_FILE = "bed";
    private static final String FREEC_DIRECTORY = "freec_dir";
    private static final String VCF_EXTENSION = "vcf_extension";
    private static final String VCF_EXTENSION_DEFAULT = ".annotation.vcf";
    private static final String CNV_RATIO_WEIGHT_FACTOR = "cnv_ratio_weight_factor";
    private static final double CNV_RATIO_WEIGHT_FACTOR_DEFAULT = 0.2;

    public static void main(final String... args) throws ParseException, IOException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);

        if (runDirectory == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Purity Ploidy Estimator (PURPLE)", options);
            System.exit(1);
        }

        final FittedCopyNumberFactory fittedCopyNumberFactory = new FittedCopyNumberFactory(
                MAX_PLOIDY,
                defaultValue(cmd, CNV_RATIO_WEIGHT_FACTOR, CNV_RATIO_WEIGHT_FACTOR_DEFAULT));

        final FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(
                MIN_PURITY,
                MAX_PURITY,
                PURITY_INCREMENTS,
                MIN_NORM_FACTOR,
                MAX_NORM_FACTOR,
                NORM_FACTOR_INCREMENTS,
                fittedCopyNumberFactory);

        LOGGER.info("Loading variant data");
        final String vcfExtention = defaultValue(cmd, VCF_EXTENSION, VCF_EXTENSION_DEFAULT);
        final VCFGermlineFile vcfFile = VCFFileLoader.loadGermlineVCF(runDirectory, vcfExtention);
        final List<GermlineVariant> variants = variants(cmd, vcfFile);
        final String refSample = vcfFile.refSample();
        final String tumorSample = vcfFile.tumorSample();

        LOGGER.info("Loading {} Freec data", tumorSample);
        final String freecDirectory = freecDirectory(cmd, runDirectory, refSample, tumorSample);
        final List<FreecCopyNumber> copyNumbers = PadCopyNumber.pad(FreecCopyNumberFactory.loadFreecCNV(freecDirectory, tumorSample));
        final List<FreecRatio> tumorRatio = FreecRatioFactory.loadTumorRatios(freecDirectory, tumorSample);
        final List<FreecRatio> normalRatio = FreecRatioFactory.loadNormalRatios(freecDirectory, tumorSample);

        LOGGER.info("Collating data");
        final EnrichedCopyNumberFactory enrichedCopyNumberFactory = new EnrichedCopyNumberFactory(
                MIN_REF_ALLELE_FREQUENCY,
                MAX_REF_ALLELE_FREQUENCY,
                MIN_COMBINED_DEPTH,
                MAX_COMBINED_DEPTH);
        final List<EnrichedCopyNumber> enrichedCopyNumbers = enrichedCopyNumberFactory.enrich(copyNumbers, variants, tumorRatio, normalRatio);

        LOGGER.info("Fitting purity");
        final List<FittedPurity> purity = fittedPurityFactory.fitPurity(enrichedCopyNumbers);
        Collections.sort(purity);

        if (!purity.isEmpty()) {
            final String cnvPath = FreecFileLoader.copyNumberPath(freecDirectory, tumorSample);
            final String purityFile = cnvPath + ".purity";
            LOGGER.info("Writing fitted purity to: {}", purityFile);
            FittedPurityWriter.writePurity(purityFile, purity);

            final String fittedFile = cnvPath + ".fitted";
            LOGGER.info("Writing fitted copy numbers to: {}", fittedFile);
            final FittedPurity bestFit = purity.get(0);
            final List<FittedCopyNumber> fittedCopyNumbers = fittedCopyNumberFactory.fittedCopyNumber(bestFit.purity(), bestFit.normFactor(), enrichedCopyNumbers);
            FittedCopyNumberWriter.writeCopyNumber(fittedFile, fittedCopyNumbers);
        }

        LOGGER.info("Complete");
    }

    private static String defaultValue(CommandLine cmd, String opt, String defaultValue) {
        return cmd.hasOption(opt) ? cmd.getOptionValue(opt) : defaultValue;
    }

    private static double defaultValue(CommandLine cmd, String opt, double defaultValue) {
        if (cmd.hasOption(opt)) {
            double result = Double.valueOf(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

    private static List<GermlineVariant> variants(CommandLine cmd, VCFGermlineFile file) throws IOException, EmptyFileException {

        final Predicate<Variant> filterPredicate = x -> x.filter().equals("PASS") || x.filter().equals(".");
        final Predicate<GenomePosition> slicerPredicate;
        if (cmd.hasOption(BED_FILE)) {
            final String bedFile = cmd.getOptionValue(BED_FILE);
            LOGGER.info("Slicing variants with bed file: " + bedFile);
            slicerPredicate = sortedSlicer(bedFile);
        } else {
            slicerPredicate = x -> true;
        }

        return file.variants().stream().filter(x -> filterPredicate.test(x) && slicerPredicate.test(x)).collect(toList());
    }

    private static String freecDirectory(CommandLine cmd, String runDirectory, String refSample, String tumorSample) {
        return cmd.hasOption(FREEC_DIRECTORY)
                ? cmd.getOptionValue(FREEC_DIRECTORY)
                : FreecFileLoader.getFreecBasePath(runDirectory, refSample, tumorSample);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, "The path containing the data for a single run");
        options.addOption(FREEC_DIRECTORY, true, "The freec data path. Defaults to ../copyNumber/sampleR_sampleT/freec/");
        options.addOption(VCF_EXTENSION, true, "VCF file extension. Defaults to " + VCF_EXTENSION_DEFAULT);
        options.addOption(BED_FILE, true, "BED file to optionally slice variants with");
        options.addOption(CNV_RATIO_WEIGHT_FACTOR, true, "CNV ratio deviation scaling");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
