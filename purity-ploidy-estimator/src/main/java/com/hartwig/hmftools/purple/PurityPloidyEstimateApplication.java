package com.hartwig.hmftools.purple;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.copynumber.freec.FreecFileLoader;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatio;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatioFactory;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatioRegions;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.FittedRegion;
import com.hartwig.hmftools.common.purple.FittedRegionFactory;
import com.hartwig.hmftools.common.purple.FittedRegionWriter;
import com.hartwig.hmftools.common.purple.ObservedRegion;
import com.hartwig.hmftools.common.purple.ObservedRegionFactory;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScoreFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurityWriter;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFGermlineFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

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

    static final double MIN_PURITY = 0.1;
    static final double MAX_PURITY = 1.0;
    static final double MIN_NORM_FACTOR = 0.33;
    static final double MAX_NORM_FACTOR = 2.0;

    private static final double MIN_REF_ALLELE_FREQUENCY = 0.4;
    private static final double MAX_REF_ALLELE_FREQUENCY = 0.65;
    private static final int MIN_COMBINED_DEPTH = 10;
    private static final int MAX_COMBINED_DEPTH = 100;
    private static final int MAX_PLOIDY = 20;
    private static final double PURITY_INCREMENTS = 0.01;
    private static final double NORM_FACTOR_INCREMENTS = 0.01;

    private static final String DEBUG = "debug";
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String FREEC_DIRECTORY = "freec_dir";
    private static final String VCF_EXTENSION = "vcf_extension";
    private static final String VCF_EXTENSION_DEFAULT = ".annotation.vcf";
    private static final String CNV_RATIO_WEIGHT_FACTOR = "cnv_ratio_weight_factor";
    private static final double CNV_RATIO_WEIGHT_FACTOR_DEFAULT = 0.2;

    public static void main(final String... args) throws ParseException, IOException, HartwigException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        final DatabaseAccess dbAccess = databaseAccess(cmd);

        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        if (runDirectory == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Purity Ploidy Estimator (PURPLE)", options);
            System.exit(1);
        }

        final FittedRegionFactory fittedRegionFactory = new FittedRegionFactory(MAX_PLOIDY,
                defaultValue(cmd, CNV_RATIO_WEIGHT_FACTOR, CNV_RATIO_WEIGHT_FACTOR_DEFAULT));

        final FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(MAX_PLOIDY, MIN_PURITY, MAX_PURITY,
                PURITY_INCREMENTS, MIN_NORM_FACTOR, MAX_NORM_FACTOR, NORM_FACTOR_INCREMENTS, fittedRegionFactory);

        LOGGER.info("Loading germline variant data");
        final String vcfExtension = defaultValue(cmd, VCF_EXTENSION, VCF_EXTENSION_DEFAULT);
        final VCFGermlineFile vcfFile = VCFFileLoader.loadGermlineVCF(runDirectory, vcfExtension);
        final List<GermlineVariant> variants = VariantFilter.passOnly(vcfFile.variants());
        final String refSample = vcfFile.refSample();
        final String tumorSample = vcfFile.tumorSample();

        LOGGER.info("Loading {} Freec data", tumorSample);
        final String freecDirectory = freecDirectory(cmd, runDirectory, refSample, tumorSample);
        final List<FreecRatio> tumorRatio = FreecRatioFactory.loadTumorRatios(freecDirectory, tumorSample);
        final List<FreecRatio> normalRatio = FreecRatioFactory.loadNormalRatios(freecDirectory, tumorSample);
        final List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(tumorRatio);

        LOGGER.info("Combining observations");
        final ObservedRegionFactory observedRegionFactory = new ObservedRegionFactory(
                MIN_REF_ALLELE_FREQUENCY, MAX_REF_ALLELE_FREQUENCY, MIN_COMBINED_DEPTH, MAX_COMBINED_DEPTH);
        final List<ObservedRegion> observedRegions = observedRegionFactory.combine(regions, variants,
                tumorRatio, normalRatio);

        LOGGER.info("Fitting purity");
        final List<FittedPurity> fittedPurities = fittedPurityFactory.fitPurity(observedRegions);
        Collections.sort(fittedPurities);

        if (!fittedPurities.isEmpty()) {
            final String purityFile = freecDirectory + File.separator + tumorSample + ".purple.purity";
            LOGGER.info("Writing fitted purity to: {}", purityFile);
            FittedPurityWriter.writePurity(purityFile, fittedPurities);

            final FittedPurity bestFit = fittedPurities.get(0);
            final List<FittedRegion> fittedRegions = fittedRegionFactory.fitRegion(bestFit.purity(),
                    bestFit.normFactor(), observedRegions);

            final List<PurpleCopyNumber> highConfidence = PurpleCopyNumberFactory.highConfidence(fittedRegions);
            final List<PurpleCopyNumber> smoothRegions = PurpleCopyNumberFactory.smooth(fittedRegions,
                    highConfidence);

            final FittedPurityScore score = FittedPurityScoreFactory.score(fittedPurities, smoothRegions);
            LOGGER.info("Persisting to database");
            dbAccess.writePurity(tumorSample, bestFit, score);
            dbAccess.writeCopynumbers(tumorSample, smoothRegions);

            if (cmd.hasOption(DEBUG)) {
                final String fittedFile = freecDirectory + File.separator + tumorSample + ".purple.fitted";
                LOGGER.info("Writing fitted copy numbers to: {}", fittedFile);
                final List<FittedRegion> broadCopyNumber = PurpleRegionZipper.insertHighConfidenceRegions(
                        highConfidence, fittedRegions);
                final List<FittedRegion> smoothCopyNumbers = PurpleRegionZipper.insertSmoothRegions(
                        smoothRegions, broadCopyNumber);
                FittedRegionWriter.writeCopyNumber(fittedFile, smoothCopyNumbers);
            }
        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static String defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt,
            @NotNull final String defaultValue) {
        return cmd.hasOption(opt) ? cmd.getOptionValue(opt) : defaultValue;
    }

    private static double defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt,
            final double defaultValue) {
        if (cmd.hasOption(opt)) {
            final double result = Double.valueOf(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

    @NotNull
    private static String freecDirectory(@NotNull final CommandLine cmd, @NotNull final String runDirectory,
            @NotNull final String refSample, @NotNull final String tumorSample) {
        return cmd.hasOption(FREEC_DIRECTORY) ?
                cmd.getOptionValue(FREEC_DIRECTORY) :
                FreecFileLoader.getFreecBasePath(runDirectory, refSample, tumorSample);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, "The path containing the data for a single run");
        options.addOption(FREEC_DIRECTORY, true,
                "The freec data path. Defaults to ../copyNumber/sampleR_sampleT/freec/");
        options.addOption(VCF_EXTENSION, true, "VCF file extension. Defaults to " + VCF_EXTENSION_DEFAULT);
        options.addOption(CNV_RATIO_WEIGHT_FACTOR, true, "CNV ratio deviation scaling");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(DEBUG, false, "Write debug fitted info to file");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static DatabaseAccess databaseAccess(CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
