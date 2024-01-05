package com.hartwig.hmftools.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleGermlineSvFile;
import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleGermlineVcfFile;
import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleSomaticSvFile;
import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleSomaticVcfFile;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.MAX_DELETED_GENES;
import static com.hartwig.hmftools.common.purple.GeneCopyNumber.listToMap;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.MemoryCalcs.calcMemoryUsage;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleSummaryData.createPurity;
import static com.hartwig.hmftools.purple.segment.Segmentation.validateObservedRegions;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MAX_SOMATIC_FIT_DELETED_PERC;
import static com.hartwig.hmftools.purple.config.PurpleConstants.TARGET_REGIONS_MAX_DELETED_GENES;
import static com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory.calculateDeletedDepthWindows;
import static com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory.validateCopyNumbers;
import static com.hartwig.hmftools.purple.fitting.BestFitFactory.buildGermlineBestFit;
import static com.hartwig.hmftools.purple.purity.FittedPurityFactory.createFittedRegionFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.AmplificationDrivers;
import com.hartwig.hmftools.common.drivercatalog.DeletionDrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.HrdData;
import com.hartwig.hmftools.common.purple.HrdDataFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.purple.fitting.PeakModelData;
import com.hartwig.hmftools.purple.fitting.SomaticPurityFitter;
import com.hartwig.hmftools.purple.gene.GeneCopyNumberBuilder;
import com.hartwig.hmftools.purple.hrd.HrdDetection;
import com.hartwig.hmftools.purple.plot.RChartData;
import com.hartwig.hmftools.purple.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.BestFit;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.purple.purity.RegionFitCalculator;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.SegmentFile;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.purple.config.AmberData;
import com.hartwig.hmftools.purple.config.CobaltData;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.config.SampleData;
import com.hartwig.hmftools.purple.config.SampleDataFiles;
import com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.purple.germline.GermlineDeletions;
import com.hartwig.hmftools.purple.germline.GermlineDrivers;
import com.hartwig.hmftools.purple.fitting.BestFitFactory;
import com.hartwig.hmftools.purple.fitting.PeakModelFile;
import com.hartwig.hmftools.purple.germline.GermlineVariants;
import com.hartwig.hmftools.purple.plot.Charts;
import com.hartwig.hmftools.purple.purity.FittedPurityFactory;
import com.hartwig.hmftools.purple.segment.Segmentation;
import com.hartwig.hmftools.purple.sv.RecoverStructuralVariants;
import com.hartwig.hmftools.purple.somatic.SomaticPeakStream;
import com.hartwig.hmftools.purple.somatic.SomaticPurityEnrichment;
import com.hartwig.hmftools.purple.somatic.SomaticStream;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;
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
        mPurpleVersion = new VersionInfo("purple.version");

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

        mSegmentation = !mConfig.DriversOnly ? new Segmentation(mReferenceData) : null;
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        try
        {
            if(mConfig.DriversOnly)
            {
                runDriversRoutine(mConfig.TumorId);
            }
            else
            {
                if(!mConfig.SampleFiles.hasValidSampleNames(mConfig))
                    System.exit(1);

                final SampleData sampleData = loadSampleData();

                if(sampleData == null)
                    System.exit(1);

                performFit(sampleData);
            }
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

        if(!mConfig.DriversOnly)
        {
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
                    new SomaticSvCache(mPurpleVersion.version(), sampleDataFiles.SomaticSvVcfFile, outputVcf, mReferenceData, mConfig)
                    : new SomaticSvCache();

            sampleData = new SampleData(referenceId, tumorId, amberData, cobaltData, svCache, somaticVariantCache);
        }
        else
        {
            sampleData = new SampleData(referenceId, tumorId, null, null, null, somaticVariantCache);
        }

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

        final CobaltChromosomes cobaltChromosomes = cobaltData.CobaltChromosomes;
        final SomaticVariantCache somaticCache = mConfig.runTumor() ? sampleData.SomaticCache : null;

        PPL_LOGGER.info("applying segmentation");
        final List<ObservedRegion> observedRegions = mSegmentation.createObservedRegions(sampleData.SvCache.variants(), amberData, cobaltData);

        if(observedRegions.isEmpty())
        {
            PPL_LOGGER.warn("no observed regions created, exiting");
            System.exit(0);
        }

        if(!validateObservedRegions(observedRegions))
        {
            PPL_LOGGER.warn("invalid observed regions, exiting");
            System.exit(0);
        }

        PPL_LOGGER.info("purple output directory: {}", mConfig.OutputDir);
        mPurpleVersion.write(mConfig.OutputDir);

        final List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        final List<ObservedRegion> fittedRegions = Lists.newArrayList();

        BestFit bestFit = null;
        PurityAdjuster purityAdjuster = null;
        Set<String> reportedGenes = Sets.newHashSet();
        SomaticStream somaticStream = null;

        RegionFitCalculator regionFitCalculator = createFittedRegionFactory(amberData.AverageTumorDepth, cobaltChromosomes, mConfig.Fitting);

        if(mConfig.runTumor())
        {
            PPL_LOGGER.info("fitting purity");

            BestFitFactory bestFitFactory = fitPurity(sampleData, observedRegions, regionFitCalculator, sampleData.SvCache.variants());

            boolean testSomaticFit = bestFitFactory.somaticFit() != null;
            bestFit = bestFitFactory.somaticFit() != null ? bestFitFactory.somaticFit() : bestFitFactory.bestNormalFit();

            purityAdjuster = new PurityAdjuster(bestFit.fit().purity(), bestFit.fit().normFactor(), cobaltChromosomes);

            buildCopyNumbers(
                    sampleData, sampleDataFiles, regionFitCalculator, purityAdjuster, observedRegions, bestFit.fit(), copyNumbers, fittedRegions);

            if(testSomaticFit)
            {
                double deletedPercent = calculateDeletedDepthWindows(copyNumbers);

                if(deletedPercent >= MAX_SOMATIC_FIT_DELETED_PERC)
                {
                    PPL_LOGGER.info(format("somatic fit(purity=%.3f ploidy=%.3f) deleted DW percent(%.3f), reverting to normal fit(purity=%.3f ploidy=%.3f)",
                            bestFit.fit().purity(), bestFit.fit().ploidy(), deletedPercent,
                            bestFitFactory.bestNormalFit().fit().purity(), bestFitFactory.bestNormalFit().fit().ploidy()));

                    // re-build using the normal fit
                    bestFit = bestFitFactory.bestNormalFit();

                    purityAdjuster = new PurityAdjuster(
                            bestFit.fit().purity(), bestFit.fit().normFactor(), cobaltChromosomes);

                    buildCopyNumbers(
                            sampleData, sampleDataFiles, regionFitCalculator, purityAdjuster, observedRegions, bestFit.fit(), copyNumbers, fittedRegions);
                }
                else
                {
                    PPL_LOGGER.debug("somatic fit deleted DW percent({})", format("%.3f", deletedPercent));
                }
            }

            sampleData.SvCache.inferMissingVariant(copyNumbers);

            geneCopyNumbers.addAll(GeneCopyNumberBuilder.createGeneCopyNumbers(
                    mReferenceData.RefGenVersion, mReferenceData.GeneTransCache, copyNumbers));

            final List<PeakModelData> somaticPeaks = Lists.newArrayList();

            PPL_LOGGER.info("modelling somatic peaks");
            final SomaticPeakStream somaticPeakStream = new SomaticPeakStream();

            final SomaticPurityEnrichment somaticPurityEnrichment = new SomaticPurityEnrichment(purityAdjuster, copyNumbers, fittedRegions);

            sampleData.SomaticCache.purityEnrich(somaticPurityEnrichment);

            List<PeakModelData> peakModelValues = somaticPeakStream.somaticPeakModel(somaticCache);
            somaticPeaks.addAll(peakModelValues);

            // at the moment the enriching of somatic variants is also contributing to the purity context, so it cannot be done afterwards
            // if the read and write process were split then so could the fitting and enriching steps
            PPL_LOGGER.info("enriching somatic variants");

            somaticStream = new SomaticStream(mConfig, mReferenceData, somaticCache, somaticPeaks);

            somaticStream.processAndWrite(purityAdjuster);

            sampleData.SvCache.write(purityAdjuster, copyNumbers, mConfig.tumorOnlyMode());

            reportedGenes.addAll(somaticStream.reportedGenes());

            FittedPurityRangeFile.write(mConfig.OutputDir, tumorId, bestFit.allFits());
            PurpleCopyNumberFile.write(PurpleCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorId), copyNumbers);
            GeneCopyNumberFile.write(GeneCopyNumberFile.generateFilenameForWriting(mConfig.OutputDir, tumorId), geneCopyNumbers);
            PeakModelFile.write(PeakModelFile.generateFilename(mConfig.OutputDir, tumorId), somaticPeaks);

            HrdData hrdData = new HrdDetection().calculateHrdData(copyNumbers);
            HrdDataFile.write(mConfig.OutputDir, tumorId, hrdData);
        }
        else
        {
            bestFit = buildGermlineBestFit();
            fittedRegions.addAll(regionFitCalculator.fitRegion(bestFit.fit().purity(), bestFit.fit().normFactor(), observedRegions));
        }

        PPL_LOGGER.info("generating QC Stats");
        final PurpleQC qcChecks = PurpleSummaryData.createQC(
                amberData.Contamination, bestFit, amberGender, cobaltGender, copyNumbers, geneCopyNumbers,
                cobaltChromosomes.germlineAberrations(), amberData.AverageTumorDepth,
                mConfig.TargetRegionsMode ? TARGET_REGIONS_MAX_DELETED_GENES : MAX_DELETED_GENES);

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
                        mPurpleVersion.version(), sampleDataFiles.GermlineSvVcfFile, mReferenceData, mConfig,
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

        if(!mConfig.germlineMode() && !mConfig.Charting.Disabled)
        {
            PPL_LOGGER.info("generating charts");

            try
            {
                Charts charts = new Charts(mConfig, mExecutorService, mReferenceData.RefGenVersion.is38());

                charts.write(
                        referenceId, tumorId, !sampleDataFiles.SomaticVcfFile.isEmpty(),
                        gender, copyNumbers, somaticStream.downsampledVariants(), sampleData.SvCache.variants(),
                        fittedRegions, Lists.newArrayList(amberData.ChromosomeBafs.values()));

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

        if(mConfig.RunDrivers)
        {
            findDrivers(tumorId, purityContext, geneCopyNumbers, somaticStream, germlineDeletions);
        }
    }

    private BestFitFactory fitPurity(
            final SampleData sampleData, final List<ObservedRegion> observedRegions,
            final RegionFitCalculator regionFitCalculator, final List<StructuralVariant> structuralVariants)
            throws ExecutionException, InterruptedException
    {
        CobaltChromosomes cobaltChromosomes = sampleData.Cobalt.CobaltChromosomes;

        List<SomaticVariant> fittingVariants = !mConfig.tumorOnlyMode() ?
                SomaticPurityFitter.findFittingVariants(sampleData.SomaticCache.variants(), observedRegions) : Lists.newArrayList();

        if(!fittingVariants.isEmpty())
        {
            PPL_LOGGER.info("somatic fitting variants({})", fittingVariants.size());
        }

        FittedPurityFactory fittedPurityFactory = new FittedPurityFactory(
                mConfig, mExecutorService, cobaltChromosomes, regionFitCalculator, observedRegions, fittingVariants);

        fittedPurityFactory.fitPurity();

        List<FittedPurity> fitCandidates = Lists.newArrayList(fittedPurityFactory.getFittedPurities());

        BestFitFactory bestFitFactory = new BestFitFactory(
                mConfig,
                sampleData.Amber.minSomaticTotalReadCount(),
                sampleData.Amber.maxSomaticTotalReadCount(),
                fitCandidates,
                fittingVariants,
                structuralVariants,
                observedRegions);

        return bestFitFactory;
    }

    private void buildCopyNumbers(
            final SampleData sampleData, final SampleDataFiles sampleDataFiles, final RegionFitCalculator regionFitCalculator,
            final PurityAdjuster purityAdjuster, final List<ObservedRegion> observedRegions, final FittedPurity fittedPurity,
            final List<PurpleCopyNumber> copyNumbers, final List<ObservedRegion> fittedRegions)
    {
        copyNumbers.clear();
        fittedRegions.clear();

        final AmberData amberData = sampleData.Amber;
        final CobaltData cobaltData = sampleData.Cobalt;

        final PurpleCopyNumberFactory copyNumberFactory = new PurpleCopyNumberFactory(
                mConfig.Fitting.MinDiploidTumorRatioCount,
                mConfig.Fitting.MinDiploidTumorRatioCountAtCentromere,
                amberData.AverageTumorDepth,
                fittedPurity.ploidy(),
                purityAdjuster,
                cobaltData.CobaltChromosomes);

        PPL_LOGGER.info("calculating copy number");
        fittedRegions.addAll(regionFitCalculator.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), observedRegions));

        copyNumberFactory.buildCopyNumbers(fittedRegions, sampleData.SvCache.variants());

        int recoveredSVCount = RecoverStructuralVariants.recoverStructuralVariants(
                sampleData, sampleDataFiles, mConfig, purityAdjuster, copyNumberFactory.copyNumbers());

        if(recoveredSVCount > 0)
        {
            PPL_LOGGER.info("reapplying segmentation with {} recovered structural variants", recoveredSVCount);
            final List<ObservedRegion> recoveredObservedRegions =
                    mSegmentation.createObservedRegions(sampleData.SvCache.variants(), amberData, cobaltData);

            PPL_LOGGER.info("recalculating copy number");
            fittedRegions.clear();
            fittedRegions.addAll(regionFitCalculator.fitRegion(fittedPurity.purity(), fittedPurity.normFactor(), recoveredObservedRegions));

            copyNumberFactory.buildCopyNumbers(fittedRegions, sampleData.SvCache.variants());
        }

        copyNumbers.addAll(copyNumberFactory.copyNumbers());

        if(!validateCopyNumbers(copyNumbers))
        {
            PPL_LOGGER.warn("invalid copy numbers, exiting");
            System.exit(0);
        }
    }

    private void runDriversRoutine(final String tumorId) throws Exception
    {
        final String purpleDataPath = mConfig.OutputDir;
        SomaticStream somaticStream = null;
        GermlineDeletions germlineDeletions = null;

        final PurityContext purityContext = PurityContextFile.read(purpleDataPath, tumorId);

        final List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(
                PurpleCopyNumberFile.generateFilenameForReading(purpleDataPath, tumorId));

        final List<ObservedRegion> fittedRegions = SegmentFile.read(SegmentFile.generateFilename(purpleDataPath, tumorId));

        final List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(
                GeneCopyNumberFile.generateFilename(purpleDataPath, tumorId));

        if(mConfig.runTumor())
        {
            String somaticVcf = purpleSomaticVcfFile(purpleDataPath, tumorId);

            SomaticVariantCache somaticVariantCache = new SomaticVariantCache(mConfig);
            ListMultimap<Chromosome, VariantHotspot> emptyHotspots = ArrayListMultimap.create(); // already annotated in VCF
            somaticVariantCache.loadSomatics(somaticVcf, emptyHotspots);

            // the counts passed in here are only for down-sampling for charting, which is not relevant for drivers
            somaticStream = new SomaticStream(mConfig, mReferenceData, somaticVariantCache, null);
            somaticStream.registerReportedVariants();
        }

        if(mConfig.runGermline())
        {
            final String germlineVcf = purpleGermlineVcfFile(purpleDataPath, tumorId);

            if(Files.exists(Paths.get(germlineVcf)))
            {
                mGermlineVariants.loadReportableVariants(germlineVcf);
                PPL_LOGGER.info("loaded {} reportable germline variants", mGermlineVariants.reportableVariants().size());
            }

            // find germline deletions since these contribute to drivers
            germlineDeletions = new GermlineDeletions(
                    mReferenceData.DriverGenes.driverGenes(), mReferenceData.GeneTransCache, mReferenceData.CohortGermlineDeletions);

            final String germlineSvVcf = purpleGermlineSvFile(purpleDataPath, tumorId);
            final List<StructuralVariant> germlineSVs = Lists.newArrayList();

            if(Files.exists(Paths.get(germlineSvVcf)))
            {
                GermlineSvCache germlineSvCache = new GermlineSvCache(
                        mPurpleVersion.version(), germlineSvVcf, mReferenceData, mConfig,
                        fittedRegions, copyNumbers, purityContext);

                germlineSVs.addAll(germlineSvCache.variants());
            }

            germlineDeletions.findDeletions(copyNumbers, fittedRegions, germlineSVs);
        }

        findDrivers(tumorId, purityContext, geneCopyNumbers, somaticStream, germlineDeletions);
    }

    private void findDrivers(
            final String tumorSample, final PurityContext purityContext,  final List<GeneCopyNumber> geneCopyNumbers,
            @Nullable final SomaticStream somaticStream, @Nullable GermlineDeletions germlineDeletions) throws IOException
    {
        final List<DriverCatalog> somaticDriverCatalog = Lists.newArrayList();
        final List<DriverCatalog> germlineDriverCatalog = Lists.newArrayList();

        PPL_LOGGER.info("generating drivers");

        final Map<String,List<GeneCopyNumber>> geneCopyNumberMap = listToMap(geneCopyNumbers);

        if(mConfig.runTumor())
        {
            List<DriverCatalog> somaticDrivers = somaticStream.buildDrivers(geneCopyNumberMap);

            somaticDriverCatalog.addAll(somaticDrivers);

            List<DriverCatalog> ampDrivers = AmplificationDrivers.findAmplifications(
                    purityContext.qc().status(), purityContext.qc().amberGender(), mReferenceData.DriverGenes,
                    purityContext.bestFit().ploidy(), geneCopyNumbers, mConfig.TargetRegionsMode);

            somaticDriverCatalog.addAll(ampDrivers);

            final DeletionDrivers delDrivers = new DeletionDrivers(purityContext.qc().status(), mReferenceData.DriverGenes);

            somaticDriverCatalog.addAll(delDrivers.deletions(geneCopyNumbers, mConfig.TargetRegionsMode));

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

    public static void main(final String... args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("Purple");

        PurpleConfig.addOptions(configBuilder);

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PurpleApplication purpleApplication = new PurpleApplication(configBuilder);
        purpleApplication.run();
    }
}
