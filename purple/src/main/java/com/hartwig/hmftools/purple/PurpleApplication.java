package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.genome.chromosome.GermlineAberration.NONE;
import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleGermlineSvFile;
import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleSomaticSvFile;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.FAIL_NO_TUMOR;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.MAX_DELETED_GENES;
import static com.hartwig.hmftools.common.purple.GeneCopyNumber.listToMap;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleSummaryData.createPurity;
import static com.hartwig.hmftools.purple.segment.Segmentation.validateObservedRegions;
import static com.hartwig.hmftools.purple.PurpleConstants.TARGET_REGIONS_MAX_DELETED_GENES;
import static com.hartwig.hmftools.purple.fitting.FittedPurityFactory.createFittedRegionFactory;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.AmplificationDrivers;
import com.hartwig.hmftools.common.driver.DeletionDrivers;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.purple.fitting.BestFit;
import com.hartwig.hmftools.purple.fitting.PurityPloidyFitter;
import com.hartwig.hmftools.purple.gene.GeneCopyNumberBuilder;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.purple.fitting.RegionFitCalculator;
import com.hartwig.hmftools.purple.germline.ChimerismDetection;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.SegmentFile;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.purple.germline.GermlineDeletions;
import com.hartwig.hmftools.purple.germline.GermlineDrivers;
import com.hartwig.hmftools.purple.fittingsnv.PeakModelFile;
import com.hartwig.hmftools.purple.germline.GermlineVariants;
import com.hartwig.hmftools.purple.plot.Charts;
import com.hartwig.hmftools.purple.segment.Segmentation;
import com.hartwig.hmftools.purple.somatic.SomaticPurityEnrichment;
import com.hartwig.hmftools.purple.somatic.SomaticStream;
import com.hartwig.hmftools.purple.somatic.SomaticVariantCache;
import com.hartwig.hmftools.purple.germline.GermlineSvCache;
import com.hartwig.hmftools.purple.sv.SomaticSvCache;

import org.jetbrains.annotations.Nullable;

public class PurpleApplication
{
    private final VersionInfo mPurpleVersion;
    private final ExecutorService mExecutorService;
    private final ReferenceData mReferenceData;
    private final PurpleConfig mConfig;

    private final GermlineVariants mGermlineVariants;
    private final Segmentation mSegmentation;

    private static final String VERSION = "version";

    private PurpleApplication(final ConfigBuilder configBuilder) throws IOException
    {
        mPurpleVersion = fromAppName(APP_NAME);

        mConfig = new PurpleConfig(mPurpleVersion.version(), configBuilder);

        if(!mConfig.isValid())
        {
            PPL_LOGGER.error("initialisation error, exiting");
            System.exit(1);
        }

        // and common reference data
        mReferenceData = new ReferenceData(configBuilder, mConfig);

        if(!mReferenceData.isValid())
        {
            PPL_LOGGER.error("invalid reference data, exiting");
            System.exit(1);
        }

        mExecutorService = Executors.newFixedThreadPool(mConfig.Threads);

        mGermlineVariants = new GermlineVariants(mConfig, mReferenceData, mPurpleVersion.version());
        mSegmentation = new Segmentation(mReferenceData);
    }

    public void run()
    {
        if(!mConfig.Fitting.hasValidValues())
        {
            PPL_LOGGER.error("invalid configured fitting values");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        try
        {
            if(!mConfig.SampleFiles.hasValidSampleNames(mConfig))
                System.exit(1);

            final SampleData sampleData = loadSampleData();

            if(sampleData == null)
                System.exit(1);

            performFit(sampleData);
        }
        catch(Throwable t)
        {
            PPL_LOGGER.error("failed processing sample({}): {}", mConfig.TumorId, t.toString());
            t.printStackTrace();
            System.exit(1);
        }

        PPL_LOGGER.info("Purple complete, mins({})", runTimeMinsStr(startTimeMs));
        mExecutorService.shutdown();
    }

    private SampleData loadSampleData() throws Exception
    {
        String referenceId = mConfig.ReferenceId;
        String tumorId = mConfig.TumorId;
        SampleDataFiles sampleDataFiles = mConfig.SampleFiles;
        SampleData sampleData = null;

        final SomaticVariantCache somaticVariantCache = new SomaticVariantCache(mConfig);

        // load amber and cobalt sample data
        final AmberData amberData = new AmberData(
                mConfig.germlineMode() ? referenceId : tumorId, sampleDataFiles.AmberDirectory, mConfig.germlineMode(),
                mReferenceData.RefGenVersion);

        final CobaltData cobaltData = new CobaltData(
                referenceId, tumorId, sampleDataFiles.CobaltDirectory, amberData.PatientGender,
                mConfig.tumorOnlyMode(), mConfig.germlineMode());

        // load structural and somatic variants
        final String outputVcf = purpleSomaticSvFile(mConfig.OutputDir, tumorId);

        final SomaticSvCache svCache = !sampleDataFiles.SomaticSvVcfFile.isEmpty() ?
                new SomaticSvCache(mPurpleVersion.version(), sampleDataFiles.SomaticSvVcfFile, outputVcf, mConfig)
                : new SomaticSvCache();

        sampleData = new SampleData(referenceId, tumorId, amberData, cobaltData, svCache, somaticVariantCache);

        if(mConfig.runTumor())
            somaticVariantCache.loadSomatics(sampleDataFiles.SomaticVcfFile, mReferenceData.SomaticHotspots);

        return sampleData;
    }

    private void performFit(final SampleData sampleData) throws Exception
    {
        String referenceId = mConfig.ReferenceId;
        String tumorId = mConfig.TumorId;
        SampleDataFiles sampleDataFiles = mConfig.SampleFiles;

        final AmberData amberData = sampleData.Amber;
        final CobaltData cobaltData = sampleData.Cobalt;

        final Gender amberGender = amberData.PatientGender;
        final Gender cobaltGender = cobaltData.gender();
        final Gender gender = mConfig.tumorOnlyMode() ? amberGender : cobaltData.gender();

        if(cobaltGender.equals(amberGender))
        {
            PPL_LOGGER.info("sample gender is {}", cobaltGender.toString().toLowerCase());
        }
        else
        {
            PPL_LOGGER.warn("Cobalt gender {} does not match Amber gender {}", cobaltGender, amberGender);
        }

        CobaltChromosomes cobaltChromosomes = cobaltData.CobaltChromosomes;
        SomaticVariantCache somaticCache = mConfig.runTumor() ? sampleData.SomaticCache : null;

        PPL_LOGGER.info("output directory: {}", mConfig.OutputDir);
        mPurpleVersion.write(mConfig.OutputDir);

        PPL_LOGGER.info("applying segmentation");
        List<ObservedRegion> observedRegions = mSegmentation.createObservedRegions(sampleData.SvCache.variants(), amberData, cobaltData);

        if(observedRegions.isEmpty() || !validateObservedRegions(observedRegions))
        {
            PPL_LOGGER.warn("no valid observed regions created, exiting");

            writeEmptyResultFiles(tumorId, sampleData, sampleDataFiles);
            System.exit(0);
        }

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        List<ObservedRegion> fittedRegions = Lists.newArrayList();

        BestFit bestFit = null;
        PurityAdjuster purityAdjuster = null;
        Set<String> reportedGenes = Sets.newHashSet();
        SomaticStream somaticStream = null;

        RegionFitCalculator regionFitCalculator = createFittedRegionFactory(amberData.AverageTumorDepth, cobaltChromosomes, mConfig.Fitting);

        if(mConfig.runTumor())
        {
            ChimerismDetection chimerismDetection = new ChimerismDetection(amberData, cobaltData, observedRegions, mReferenceData.RefGenVersion);
            chimerismDetection.run();

            PPL_LOGGER.info("fitting purity");

            PurityPloidyFitter purityPloidyFitter = new PurityPloidyFitter(
                    mConfig, mReferenceData, sampleData, mExecutorService, regionFitCalculator, observedRegions, gender);

            purityPloidyFitter.run();

            fittedRegions.addAll(purityPloidyFitter.fittedRegions());

            bestFit = purityPloidyFitter.finalFit();

            purityAdjuster = purityPloidyFitter.purityAdjuster();

            copyNumbers.addAll(purityPloidyFitter.copyNumbers());

            sampleData.SvCache.inferMissingVariant(copyNumbers);

            geneCopyNumbers.addAll(GeneCopyNumberBuilder.createGeneCopyNumbers(
                    mReferenceData.RefGenVersion, mReferenceData.GeneTransCache, copyNumbers));

            SomaticPurityEnrichment somaticPurityEnrichment = new SomaticPurityEnrichment(purityAdjuster, copyNumbers, fittedRegions);

            sampleData.SomaticCache.purityEnrich(somaticPurityEnrichment);

            // at the moment the enriching of somatic variants is also contributing to the purity context, so it cannot be done afterwards
            // if the read and write process were split then so could the fitting and enriching steps
            PPL_LOGGER.info("enriching somatic variants");

            somaticStream = new SomaticStream(mConfig, mReferenceData, somaticCache);

            somaticStream.processAndWrite(purityAdjuster);

            sampleData.SvCache.write(purityAdjuster, copyNumbers, mConfig.tumorOnlyMode(), amberGender);

            reportedGenes.addAll(somaticStream.reportedGenes());

            FittedPurityRangeFile.write(mConfig.OutputDir, tumorId, bestFit.AllFits);
            PurpleCopyNumberFile.write(PurpleCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorId), copyNumbers);
            GeneCopyNumberFile.write(GeneCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorId), geneCopyNumbers);
            PeakModelFile.write(PeakModelFile.generateFilename(mConfig.OutputDir, tumorId), somaticStream.peakModelData());
        }
        else
        {
            bestFit = PurityPloidyFitter.buildGermlineBestFit();
            fittedRegions.addAll(regionFitCalculator.fitRegion(bestFit.Fit.purity(), bestFit.Fit.normFactor(), observedRegions));
        }

        PPL_LOGGER.debug("generating QC stats");

        final PurpleQC qcChecks = PurpleSummaryData.createQC(
                amberData.Contamination, bestFit, amberGender, cobaltGender, copyNumbers, geneCopyNumbers,
                cobaltChromosomes.germlineAberrations(), amberData.AverageTumorDepth,
                mConfig.TargetRegionsMode ? TARGET_REGIONS_MAX_DELETED_GENES : MAX_DELETED_GENES,
                somaticCache != null ? somaticCache.tincLevel() : 0);

        final PurityContext purityContext = createPurity(bestFit, gender, mConfig, qcChecks, copyNumbers, somaticStream, sampleData.SvCache);

        PurityContextFile.write(mConfig.OutputDir, tumorId, purityContext);
        SegmentFile.write(SegmentFile.generateFilename(mConfig.OutputDir, tumorId), fittedRegions);

        GermlineDeletions germlineDeletions = null;

        if(mConfig.runGermline())
        {
            mGermlineVariants.processAndWrite(
                    referenceId, tumorId, sampleDataFiles.GermlineVcfFile, purityAdjuster, copyNumbers, reportedGenes);

            GermlineSvCache germlineSvCache;

            if(!sampleDataFiles.GermlineSvVcfFile.isEmpty())
            {
                germlineSvCache = new GermlineSvCache(
                        mPurpleVersion.version(), sampleDataFiles.GermlineSvVcfFile, mConfig,
                        fittedRegions, copyNumbers, purityContext);

                germlineSvCache.write(purpleGermlineSvFile(mConfig.OutputDir, tumorId));
            }
            else
            {
                germlineSvCache = new GermlineSvCache();
            }

            germlineDeletions = new GermlineDeletions(
                    mReferenceData.DriverGenes.driverGenes(), mReferenceData.GeneTransCache, mReferenceData.CohortGermlineDeletions);

            germlineDeletions.findDeletions(copyNumbers, fittedRegions, germlineSvCache.variants());

            GermlineDeletion.write(GermlineDeletion.generateFilename(mConfig.OutputDir, tumorId), germlineDeletions.getDeletions());
        }

        List<DriverSourceData> driverSourceData = Lists.newArrayList();

        if(mConfig.RunDrivers)
        {
            findDrivers(tumorId, purityContext, geneCopyNumbers, somaticStream, germlineDeletions, driverSourceData);
        }

        if(!mConfig.germlineMode() && !mConfig.Charting.Disabled)
        {
            PPL_LOGGER.info("generating charts");

            try
            {
                Charts charts = new Charts(mConfig, mExecutorService, mReferenceData.RefGenVersion.is38());

                List<ObservedRegion> regionsForVis = fittedRegions.stream()
                        .filter(x -> x.germlineStatus() != GermlineStatus.EXCLUDED).collect(Collectors.toList());

                charts.write(
                        referenceId, tumorId, !sampleDataFiles.SomaticVcfFile.isEmpty(),
                        gender, copyNumbers, somaticStream.plottingVariants(), sampleData.SvCache.variants(),
                        regionsForVis, Lists.newArrayList(amberData.ChromosomeBafs.values()), driverSourceData);
                
                // clean up any temporary files
                // RChartData.cleanupFiles(mConfig, tumorId);
            }
            catch(Exception e)
            {
                PPL_LOGGER.error("charting error: {}", e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    private void findDrivers(
            final String tumorSample, final PurityContext purityContext,  final List<GeneCopyNumber> geneCopyNumbers,
            @Nullable final SomaticStream somaticStream, @Nullable GermlineDeletions germlineDeletions,
            final List<DriverSourceData> driverSourceData) throws IOException
    {
        List<DriverCatalog> somaticDriverCatalog = Lists.newArrayList();
        List<DriverCatalog> germlineDriverCatalog = Lists.newArrayList();

        PPL_LOGGER.info("generating drivers");

        Map<String,List<GeneCopyNumber>> geneCopyNumberMap = listToMap(geneCopyNumbers);

        if(mConfig.runTumor())
        {
            List<DriverCatalog> somaticDrivers = somaticStream.buildDrivers(geneCopyNumberMap, driverSourceData);

            somaticDriverCatalog.addAll(somaticDrivers);

            List<DriverCatalog> ampDrivers = AmplificationDrivers.findAmplifications(
                    purityContext.qc().status(), purityContext.qc().amberGender(), mReferenceData.DriverGenes,
                    purityContext.bestFit().ploidy(), geneCopyNumbers, mConfig.TargetRegionsMode);

            addGeneDriverSourceData(ampDrivers, geneCopyNumbers, driverSourceData);

            somaticDriverCatalog.addAll(ampDrivers);

            DeletionDrivers delDriverFinder = new DeletionDrivers(purityContext.qc().status(), mReferenceData.DriverGenes);

            List<DriverCatalog> delDrivers = delDriverFinder.deletions(geneCopyNumbers, mConfig.TargetRegionsMode);

            somaticDriverCatalog.addAll(delDrivers);

            addGeneDriverSourceData(delDrivers, geneCopyNumbers, driverSourceData);

            DriverCatalogFile.write(DriverCatalogFile.generateFilenameForWriting(mConfig.OutputDir, tumorSample, true), somaticDriverCatalog);
        }

        if(mConfig.runGermline())
        {
            final GermlineDrivers germlineDrivers = new GermlineDrivers(mReferenceData.DriverGenes.driverGenes());
            germlineDriverCatalog.addAll(germlineDrivers.findDrivers(mGermlineVariants.reportableVariants(), geneCopyNumberMap));

            if(germlineDeletions != null)
                germlineDriverCatalog.addAll(germlineDeletions.getDrivers());

            DriverCatalogFile.write(DriverCatalogFile.generateFilenameForWriting(mConfig.OutputDir, tumorSample, false), germlineDriverCatalog);
        }
    }

    private void addGeneDriverSourceData(
            final List<DriverCatalog> geneDivers, final List<GeneCopyNumber> geneCopyNumbers, final List<DriverSourceData> driverSourceData)
    {
        for(DriverCatalog driverRecord : geneDivers)
        {
            GeneCopyNumber geneCopyNumber = geneCopyNumbers.stream()
                    .filter(x -> x.geneName().equals(driverRecord.gene())).findFirst().orElse(null);

            if(geneCopyNumber != null)
            {
                driverSourceData.add(new DriverSourceData(driverRecord, geneCopyNumber));
            }
        }
    }

    private void writeEmptyResultFiles(
            final String tumorId, final SampleData sampleData, final SampleDataFiles sampleDataFiles) throws IOException
    {
        Gender gender = Gender.FEMALE;

        if(mConfig.runTumor())
        {
            SomaticStream somaticStream = new SomaticStream(mConfig, mReferenceData, sampleData.SomaticCache);
            somaticStream.processAndWrite(null);

            sampleData.SvCache.write(null, Collections.emptyList(), mConfig.tumorOnlyMode(), gender);

            FittedPurityRangeFile.write(mConfig.OutputDir, tumorId, Collections.emptyList());
            PurpleCopyNumberFile.write(PurpleCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorId), Collections.emptyList());
            GeneCopyNumberFile.write(GeneCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorId), Collections.emptyList());
            PeakModelFile.write(PeakModelFile.generateFilename(mConfig.OutputDir, tumorId), Collections.emptyList());

            FittedPurity fittedPurity = ImmutableFittedPurity.builder()
                    .purity(0).ploidy(0).score(0).diploidProportion(0).normFactor(0).somaticPenalty(0).build();

            FittedPurityScore fittedPurityScore = ImmutableFittedPurityScore.builder()
                    .maxPurity(0).minPurity(0).maxPloidy(0).minPloidy(0).maxDiploidProportion(0).minDiploidProportion(0).build();

            PurpleQC purpleQC = ImmutablePurpleQC.builder()
                    .method(FittedPurityMethod.NO_TUMOR).purity(0).contamination(0).cobaltGender(gender)
                    .unsupportedCopyNumberSegments(0).deletedGenes(0).amberGender(gender).lohPercent(0).copyNumberSegments(0)
                    .status(List.of(FAIL_NO_TUMOR)).germlineAberrations(List.of(NONE)).amberMeanDepth(0).tincLevel(0).build();

            PurityContext purityContext =  ImmutablePurityContext.builder()
                    .bestFit(fittedPurity)
                    .method(FittedPurityMethod.NO_TUMOR)
                    .gender(gender)
                    .runMode(mConfig.runMode())
                    .score(fittedPurityScore)
                    .polyClonalProportion(0)
                    .wholeGenomeDuplication(false)
                    .microsatelliteIndelsPerMb(0)
                    .microsatelliteStatus(MicrosatelliteStatus.UNKNOWN)
                    .tumorMutationalLoad(0)
                    .tumorMutationalLoadStatus(TumorMutationalStatus.UNKNOWN)
                    .tumorMutationalBurdenPerMb(0)
                    .tumorMutationalBurdenStatus(TumorMutationalStatus.UNKNOWN)
                    .svTumorMutationalBurden(0)
                    .qc(purpleQC)
                    .targeted(mConfig.TargetRegionsMode)
                    .build();

            PurityContextFile.write(mConfig.OutputDir, tumorId, purityContext);
            SegmentFile.write(SegmentFile.generateFilename(mConfig.OutputDir, tumorId), Collections.emptyList());

            DriverCatalogFile.write(DriverCatalogFile.generateFilenameForWriting(
                    mConfig.OutputDir, tumorId, true), Collections.emptyList());
        }

        if(mConfig.runGermline())
        {
            if(!sampleDataFiles.GermlineSvVcfFile.isEmpty())
            {
                GermlineSvCache germlineSvCache = new GermlineSvCache(
                        mPurpleVersion.version(), sampleDataFiles.GermlineSvVcfFile, mConfig,
                        Collections.emptyList(), Collections.emptyList(), null);

                germlineSvCache.write(purpleGermlineSvFile(mConfig.OutputDir, tumorId));
            }

            GermlineDeletion.write(GermlineDeletion.generateFilename(mConfig.OutputDir, tumorId), Collections.emptyList());
            DriverCatalogFile.write(DriverCatalogFile.generateFilenameForWriting(mConfig.OutputDir, tumorId, false), Collections.emptyList());
        }
    }

    private static final String APP_NAME = "Purple";

    public static void main(final String... args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        PurpleConfig.registerConfig(configBuilder);

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PurpleApplication purpleApplication = new PurpleApplication(configBuilder);
        purpleApplication.run();
    }
}
