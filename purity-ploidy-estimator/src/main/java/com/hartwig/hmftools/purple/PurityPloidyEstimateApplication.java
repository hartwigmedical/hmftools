package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purple.purity.FittedPurityScoreFactory.polyclonalProportion;
import static com.hartwig.hmftools.patientdb.LoadPurpleData.persistToDatabase;
import static com.hartwig.hmftools.purple.PurpleRegionZipper.updateRegionsWithCopyNumbers;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFactory;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.BestFit;
import com.hartwig.hmftools.common.purple.purity.BestFitFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFactory;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.purity.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFactory;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactoryV2;
import com.hartwig.hmftools.common.purple.region.FittedRegionFile;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.filter.SGTFilter;
import com.hartwig.hmftools.common.variant.recovery.RecoverStructuralVariants;
import com.hartwig.hmftools.common.version.VersionInfo;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.purple.config.CircosConfig;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;
import com.hartwig.hmftools.purple.config.DBConfig;
import com.hartwig.hmftools.purple.config.FitScoreConfig;
import com.hartwig.hmftools.purple.config.FittingConfig;
import com.hartwig.hmftools.purple.config.SmoothingConfig;
import com.hartwig.hmftools.purple.config.SomaticConfig;
import com.hartwig.hmftools.purple.config.StructuralVariantConfig;
import com.hartwig.hmftools.purple.plot.ChartWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class PurityPloidyEstimateApplication {

    private static final Logger LOGGER = LogManager.getLogger(PurityPloidyEstimateApplication.class);

    private static final int THREADS_DEFAULT = 2;

    private static final String THREADS = "threads";
    private static final String VERSION = "version";
    private static final String EXPERIMENTAL = "experimental";

    public static void main(final String... args)
            throws ParseException, IOException, SQLException, ExecutionException, InterruptedException {
        new PurityPloidyEstimateApplication(args);
    }

    private PurityPloidyEstimateApplication(final String... args)
            throws ParseException, IOException, SQLException, ExecutionException, InterruptedException {
        final VersionInfo version = new VersionInfo("purple.version");
        LOGGER.info("PURPLE version: {}", version.version());

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (cmd.hasOption(VERSION)) {
            System.exit(0);
        }

        final int threads = cmd.hasOption(THREADS) ? Integer.valueOf(cmd.getOptionValue(THREADS)) : THREADS_DEFAULT;
        final ExecutorService executorService = Executors.newFixedThreadPool(threads);
        try {
            // JOBA: Get common config
            final ConfigSupplier configSupplier = new ConfigSupplier(cmd, options);
            final CommonConfig config = configSupplier.commonConfig();
            final String outputDirectory = config.outputDirectory();
            final String tumorSample = config.tumorSample();

            // JOBA: Read Gene Panel
            final List<HmfTranscriptRegion> genePanel = HmfGenePanelSupplier.allGeneList();

            // JOBA: Load Amber Data
            final Gender amberGender = configSupplier.amberData().gender();
            final Multimap<Chromosome, AmberBAF> bafs = configSupplier.amberData().bafs();
            int averageTumorDepth = configSupplier.amberData().averageTumorDepth();

            // JOBA: Load Cobalt Data
            final Gender cobaltGender = configSupplier.cobaltData().gender();
            final ListMultimap<Chromosome, CobaltRatio> ratios = configSupplier.cobaltData().ratios();

            // JOBA: Gender
            if (cobaltGender.equals(amberGender)) {
                LOGGER.info("Sample gender is {}", cobaltGender.toString().toLowerCase());
            } else {
                LOGGER.warn("COBALT gender {} does not match AMBER gender {}", cobaltGender, amberGender);
            }

            // JOBA: Load structural and somatic variants
            final PurpleStructuralVariantSupplier structuralVariants = structuralVariants(configSupplier);
            final List<SomaticVariant> allSomatics = somaticVariants(configSupplier);
            final List<SomaticVariant> snpSomatics = allSomatics.stream().filter(SomaticVariant::isSnp).collect(Collectors.toList());

            LOGGER.info("Applying segmentation");
            final Segmentation segmentation = new Segmentation(configSupplier, cobaltGender, ratios, bafs);
            final List<ObservedRegion> observedRegions = segmentation.createSegments(structuralVariants.variants());

            LOGGER.info("Fitting purity");
            final FitScoreConfig fitScoreConfig = configSupplier.fitScoreConfig();
            final FittedRegionFactory fittedRegionFactory = createFittedRegionFactory(averageTumorDepth, cobaltGender, fitScoreConfig);
            final BestFit bestFit =
                    fitPurity(executorService, configSupplier, cobaltGender, snpSomatics, observedRegions, fittedRegionFactory);
            final FittedPurity fittedPurity = bestFit.fit();
            final PurityAdjuster purityAdjuster = new PurityAdjuster(cobaltGender, fittedPurity);

            final SmoothingConfig smoothingConfig = configSupplier.smoothingConfig();
            final PurpleCopyNumberFactory copyNumberFactory = new PurpleCopyNumberFactory(smoothingConfig.minDiploidTumorRatioCount(),
                    smoothingConfig.minDiploidTumorRatioCountAtCentromere(),
                    purityAdjuster);

            LOGGER.info("Calculating copy number");
            List<FittedRegion> fittedRegions =
                    fittedRegionFactory.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), observedRegions);
            copyNumberFactory.invoke(fittedRegions, structuralVariants.variants());

            StructuralVariantConfig svConfig = configSupplier.structuralVariantConfig();
            if (svConfig.recoveryFile().isPresent()) {
                final String vcfRecoveryFile = svConfig.recoveryFile().get().toString();
                try (final RecoverStructuralVariants recovery = new RecoverStructuralVariants(purityAdjuster,
                        vcfRecoveryFile,
                        copyNumberFactory.copyNumbers())) {
                    final Collection<VariantContext> recoveredVariants = recovery.recoverVariants(structuralVariants.variants());
                    if (!recoveredVariants.isEmpty()) {
                        recoveredVariants.forEach(structuralVariants::recoverVariant);

                        LOGGER.info("Recalculating segmentation with {} recovered structural variants", recoveredVariants.size());
                        final List<ObservedRegion> recoveredObservedRegions = segmentation.createSegments(structuralVariants.variants());

                        LOGGER.info("Recalculating copy number");
                        fittedRegions = fittedRegionFactory.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), recoveredObservedRegions);
                        copyNumberFactory.invoke(fittedRegions, structuralVariants.variants());
                    }
                }
            }

            final List<PurpleCopyNumber> copyNumbers = copyNumberFactory.copyNumbers();
            final List<PurpleCopyNumber> germlineDeletions = copyNumberFactory.germlineDeletions();
            final List<FittedRegion> enrichedFittedRegions = updateRegionsWithCopyNumbers(fittedRegions, copyNumbers);

            final PurityContext purityContext = ImmutablePurityContext.builder()
                    .version(version.version())
                    .bestFit(bestFit.fit())
                    .status(bestFit.status())
                    .gender(cobaltGender)
                    .score(bestFit.score())
                    .polyClonalProportion(polyclonalProportion(copyNumbers))
                    .build();

            final List<PurityAdjustedSomaticVariant> enrichedSomatics =
                    new PurityAdjustedSomaticVariantFactory(purityAdjuster, copyNumbers, enrichedFittedRegions).create(allSomatics);

            final List<GeneCopyNumber> geneCopyNumbers =
                    GeneCopyNumberFactory.geneCopyNumbers(genePanel, copyNumbers, germlineDeletions, enrichedSomatics);

            LOGGER.info("Generating QC Stats");
            final PurpleQC qcChecks = PurpleQCFactory.create(bestFit.fit(), copyNumbers, amberGender, cobaltGender, geneCopyNumbers);

            final DBConfig dbConfig = configSupplier.dbConfig();
            if (dbConfig.enabled()) {
                final DatabaseAccess dbAccess = databaseAccess(dbConfig);
                persistToDatabase(dbAccess,
                        tumorSample,
                        bestFit.bestFitPerPurity(),
                        copyNumbers,
                        germlineDeletions,
                        enrichedFittedRegions,
                        purityContext,
                        qcChecks,
                        geneCopyNumbers);
            }

            LOGGER.info("Writing purple data to: {}", outputDirectory);
            version.write(outputDirectory);
            PurpleQCFile.write(PurpleQCFile.generateFilename(outputDirectory, tumorSample), qcChecks);
            FittedPurityFile.write(outputDirectory, tumorSample, purityContext);
            FittedPurityRangeFile.write(outputDirectory, tumorSample, bestFit.bestFitPerPurity());
            PurpleCopyNumberFile.write(PurpleCopyNumberFile.generateFilename(outputDirectory, tumorSample), copyNumbers);
            PurpleCopyNumberFile.write(PurpleCopyNumberFile.generateGermlineFilename(outputDirectory, tumorSample), germlineDeletions);
            FittedRegionFile.write(FittedRegionFile.generateFilename(outputDirectory, tumorSample), enrichedFittedRegions);
            GeneCopyNumberFile.write(GeneCopyNumberFile.generateFilename(outputDirectory, tumorSample), geneCopyNumbers);
            structuralVariants.write();

            final CircosConfig circosConfig = configSupplier.circosConfig();
            LOGGER.info("Writing plots to: {}", circosConfig.plotDirectory());
            new ChartWriter(tumorSample, circosConfig.plotDirectory()).write(purityContext.bestFit(),
                    purityContext.score(),
                    copyNumbers,
                    enrichedSomatics);

            LOGGER.info("Writing circos data to: {}", circosConfig.circosDirectory());
            new GenerateCircosData(configSupplier, executorService).write(cobaltGender,
                    copyNumbers,
                    enrichedSomatics,
                    structuralVariants.variants(),
                    fittedRegions,
                    Lists.newArrayList(bafs.values()));
        } finally {
            executorService.shutdown();
        }
        LOGGER.info("Complete");
    }

    private BestFit fitPurity(final ExecutorService executorService, final ConfigSupplier configSupplier, final Gender cobaltGender,
            final List<SomaticVariant> snpSomatics, final List<ObservedRegion> observedRegions,
            final FittedRegionFactory fittedRegionFactory) throws ExecutionException, InterruptedException {
        final FittingConfig fittingConfig = configSupplier.fittingConfig();
        final SomaticConfig somaticConfig = configSupplier.somaticConfig();
        final FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(executorService,
                cobaltGender,
                fittingConfig.maxPloidy(),
                fittingConfig.minPurity(),
                fittingConfig.maxPurity(),
                fittingConfig.purityIncrement(),
                fittingConfig.minNormFactor(),
                fittingConfig.maxNormFactor(),
                fittingConfig.normFactorIncrement(),
                somaticConfig.somaticDeviationWeight(),
                fittedRegionFactory,
                observedRegions,
                snpSomatics);

        final List<FittedPurity> bestFitPerPurity = fittedPurityFactory.bestFitPerPurity();

        final BestFitFactory bestFitFactory = new BestFitFactory(somaticConfig.minTotalVariants(),
                somaticConfig.minPeakVariants(),
                somaticConfig.highlyDiploidPercentage(),
                somaticConfig.minSomaticPurity(),
                somaticConfig.minSomaticPuritySpread(),
                bestFitPerPurity,
                snpSomatics);
        return bestFitFactory.bestFit();
    }

    @NotNull
    private FittedRegionFactory createFittedRegionFactory(final int averageTumorDepth, final Gender cobaltGender,
            final FitScoreConfig fitScoreConfig) {
        return new FittedRegionFactoryV2(cobaltGender,
                averageTumorDepth,
                fitScoreConfig.ploidyPenaltyFactor(),
                fitScoreConfig.ploidyPenaltyStandardDeviation(),
                fitScoreConfig.ploidyPenaltyMinStandardDeviationPerPloidy(),
                fitScoreConfig.ploidyPenaltyMajorAlleleSubOneMultiplier(),
                fitScoreConfig.ploidyPenaltyMajorAlleleSubOneAdditional(),
                fitScoreConfig.ploidyPenaltyBaselineDeviation());
    }

    @NotNull
    private static PurpleStructuralVariantSupplier structuralVariants(@NotNull final ConfigSupplier configSupplier) throws IOException {
        final CommonConfig commonConfig = configSupplier.commonConfig();
        final StructuralVariantConfig svConfig = configSupplier.structuralVariantConfig();
        if (svConfig.file().isPresent()) {
            final String filePath = svConfig.file().get().toString();
            final String outputPath = commonConfig.outputDirectory() + File.separator + commonConfig.tumorSample() + ".purple.sv.vcf.gz";

            LOGGER.info("Loading structural variants from {}", filePath);
            return new PurpleStructuralVariantSupplier(filePath, outputPath);
        } else {
            return new PurpleStructuralVariantSupplier();
        }
    }

    @NotNull
    private static List<SomaticVariant> somaticVariants(@NotNull final ConfigSupplier configSupplier) throws IOException {
        final SomaticConfig config = configSupplier.somaticConfig();
        if (config.file().isPresent()) {
            String filename = config.file().get().toString();
            LOGGER.info("Loading somatic variants from {}", filename);

            SomaticVariantFactory factory =
                    SomaticVariantFactory.filteredInstance(new PassingVariantFilter(), new NTFilter(), new SGTFilter());

            return factory.fromVCFFile(configSupplier.commonConfig().tumorSample(), filename);
        } else {
            LOGGER.info("Somatic variants support disabled.");
            return Collections.emptyList();
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        ConfigSupplier.addOptions(options);

        options.addOption(THREADS, true, "Number of threads (default 2)");
        options.addOption(EXPERIMENTAL, false, "Anything goes!");
        options.addOption(VERSION, false, "Exit after displaying version info.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final DBConfig dbConfig) throws SQLException {
        return new DatabaseAccess(dbConfig.user(), dbConfig.password(), dbConfig.url());
    }
}
