package com.hartwig.hmftools.orange.algo;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxData;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxDataLoader;
import com.hartwig.hmftools.orange.algo.linx.LinxData;
import com.hartwig.hmftools.orange.algo.linx.LinxDataLoader;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.flagstat.Flagstat;
import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.metrics.WGSMetrics;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangePlots;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeSample;
import com.hartwig.hmftools.datamodel.orange.OrangePlots;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeSample;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.wildtype.WildTypeGene;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeRnaConfig;
import com.hartwig.hmftools.orange.OrangeWGSRefConfig;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaDataFactory;
import com.hartwig.hmftools.orange.algo.immuno.ImmuneEscapeInterpreter;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpreter;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpreter;
import com.hartwig.hmftools.orange.algo.linx.LinxReportableClusters;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.plot.DummyPlotManager;
import com.hartwig.hmftools.orange.algo.plot.FileBasedPlotManager;
import com.hartwig.hmftools.orange.algo.plot.PlotManager;
import com.hartwig.hmftools.orange.algo.purple.GermlineGainDeletionFactory;
import com.hartwig.hmftools.orange.algo.purple.PurpleData;
import com.hartwig.hmftools.orange.algo.purple.PurpleDataLoader;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpreter;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariantFactory;
import com.hartwig.hmftools.orange.algo.sage.GermlineMVLHFactory;
import com.hartwig.hmftools.orange.algo.sigs.SigsEtiologiesLoader;
import com.hartwig.hmftools.orange.algo.sigs.SigsInterpreter;
import com.hartwig.hmftools.orange.algo.util.ReportLimiter;
import com.hartwig.hmftools.orange.algo.virus.VirusInterpreter;
import com.hartwig.hmftools.orange.algo.wildtype.WildTypeAlgo;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableSample;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapper;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapping;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMappingFile;
import com.hartwig.hmftools.orange.cohort.mapping.DoidCohortMapper;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.OrangeConversion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class OrangeAlgo
{
    private final DoidEntry mDoidEntry;
    private final CohortMapper mCohortMapper;
    private final Map<String,DriverGene> mDriverGenes;
    private final Map<String, String> mEtiologyPerSignature;
    private final EnsemblDataCache mEnsemblDataCache;
    private final PlotManager mPlotManager;

    private boolean mSuppressGeneWarnings;

    public static OrangeAlgo fromConfig(final OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading DOID database from {}", config.doidJsonFile());
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJsonFile());
        DoidParents doidParentModel = DoidParents.fromEdges(doidEntry.edges());

        LOGGER.info("Reading cohort mappings from {}", config.cohortMappingTsv());
        List<CohortMapping> mappings = CohortMappingFile.read(config.cohortMappingTsv());
        LOGGER.info(" Reading {} cohort mappings", mappings.size());
        CohortMapper mapper = new DoidCohortMapper(doidParentModel, mappings);

        LOGGER.info("Reading driver genes from {}", config.driverGenePanelTsv());
        List<DriverGene> driverGenes = DriverGeneFile.read(config.driverGenePanelTsv());
        LOGGER.info(" Read {} driver genes", driverGenes.size());

        LOGGER.info("Reading signatures etiology from {}", config.signaturesEtiologyTsv());
        Map<String, String> etiologyPerSignature = SigsEtiologiesLoader.read(config.signaturesEtiologyTsv());
        LOGGER.info(" Read {} signatures etiology", etiologyPerSignature.size());

        LOGGER.info("Reading ensembl data cache from {}", config.ensemblDataDirectory());
        EnsemblDataCache ensemblDataCache = loadEnsemblDataCache(config);
        LOGGER.info(" Read ensembl data dir");

        String outputDir = config.outputDir();
        PlotManager plotManager = !outputDir.isEmpty() ? new FileBasedPlotManager(outputDir) : new DummyPlotManager();

        return new OrangeAlgo(doidEntry, mapper, driverGenes, etiologyPerSignature, ensemblDataCache, plotManager);
    }

    private OrangeAlgo(
            final DoidEntry doidEntry, final CohortMapper cohortMapper, final List<DriverGene> driverGenes,
            final Map<String, String> etiologyPerSignature, final EnsemblDataCache ensemblDataCache, final PlotManager plotManager)
    {
        mDoidEntry = doidEntry;
        mCohortMapper = cohortMapper;

        mDriverGenes = Maps.newHashMap();
        driverGenes.forEach(x -> mDriverGenes.put(x.gene(), x));

        mEtiologyPerSignature = etiologyPerSignature;
        mEnsemblDataCache = ensemblDataCache;
        mPlotManager = plotManager;
        mSuppressGeneWarnings = false;
    }

    public OrangeRecord run(final OrangeConfig config) throws Exception
    {
        Set<DoidNode> configuredPrimaryTumor = loadConfiguredPrimaryTumor(config);
        String pipelineVersion = determinePipelineVersion(config);
        OrangeSample refSample = loadSampleData(config, false);
        OrangeSample tumorSample = loadSampleData(config, true);

        PurpleData purpleData = loadPurpleData(config);
        LinxData linxData = loadLinxData(config);
        Map<String, Double> mvlhPerGene = loadGermlineMVLHPerGene(config, mDriverGenes);
        ChordData chord = loadChordAnalysis(config);
        LilacSummaryData lilac = loadLilacData(config);
        VirusInterpreterData virusInterpreter = loadVirusInterpreterData(config);
        CuppaData cuppa = loadCuppaData(config);
        List<PeachGenotype> peach = loadPeachData(config);
        List<SignatureAllocation> sigAllocations = loadSigAllocations(config);
        IsofoxData isofoxData = loadIsofoxData(config);

        LinxInterpreter linxInterpreter = new LinxInterpreter(mEnsemblDataCache);

        LinxRecord linx = linxInterpreter.interpret(linxData);

        PaveAlgo pave = new PaveAlgo(mEnsemblDataCache, !mSuppressGeneWarnings);

        PurpleVariantFactory purpleVariantFactory = new PurpleVariantFactory(pave);

        GermlineGainDeletionFactory germlineGainDeletionFactory = new GermlineGainDeletionFactory(mEnsemblDataCache);

        List<DriverGene> driverGenes = new ArrayList<>(mDriverGenes.values());

        PurpleInterpreter purpleInterpreter = new PurpleInterpreter(purpleVariantFactory, germlineGainDeletionFactory);

        PurpleRecord purple = purpleInterpreter.interpret(purpleData);

        ImmuneEscapeRecord immuneEscape = ImmuneEscapeInterpreter.interpret(purple, linx);

        IsofoxRecord isofox = null;
        if(isofoxData != null)
        {
            IsofoxInterpreter isofoxInterpreter = new IsofoxInterpreter(driverGenes, linx);
            isofox = isofoxInterpreter.interpret(isofoxData);
        }

        List<WildTypeGene> wildTypeGenes = Lists.newArrayList();

        if(WildTypeAlgo.wildTypeCallingAllowed(purple.fit().qc().status()))
        {
            wildTypeGenes = WildTypeAlgo.determineWildTypeGenes(
                    mDriverGenes,
                    purple.somaticVariants(), purple.germlineVariants(), purple.somaticGainsDels(),
                    linx.fusions(), linx.somaticHomozygousDisruptions(), linx.somaticBreakends());

            LOGGER.info("Identified {} of {} driver genes to be wild-type", wildTypeGenes.size(), driverGenes.size());
        }
        else
        {
            LOGGER.info("Wild-type calling skipped due to insufficient tumor sample quality");
        }

        boolean hasRefSample = config.wgsRefConfig() != null && config.wgsRefConfig().referenceSampleId() != null;

        OrangeRecord report = ImmutableOrangeRecord.builder()
                .sampleId(config.tumorSampleId())
                .samplingDate(config.samplingDate())
                .experimentType(config.experimentType())
                .configuredPrimaryTumor(ConversionUtil.mapToIterable(configuredPrimaryTumor, OrangeConversion::convert))
                .refGenomeVersion(config.refGenomeVersion())
                .pipelineVersion(pipelineVersion)
                .refSample(refSample)
                .tumorSample(tumorSample)
                .germlineMVLHPerGene(mvlhPerGene)
                .purple(purple)
                .linx(linx)
                .wildTypeGenes(wildTypeGenes)
                .isofox(isofox)
                .lilac(lilac != null ? OrangeConversion.convert(lilac, hasRefSample, config.rnaConfig() != null) : null)
                .immuneEscape(immuneEscape)
                .virusInterpreter(virusInterpreter != null ? VirusInterpreter.interpret(virusInterpreter) : null)
                .chord(chord != null ? OrangeConversion.convert(chord) : null)
                .cuppa(cuppa)
                .peach(ConversionUtil.mapToIterable(peach, OrangeConversion::convert))
                .sigAllocations(SigsInterpreter.interpret(sigAllocations, mEtiologyPerSignature))
                .plots(buildPlots(config))
                .build();

        verifyPlots(report.plots(), linxData);

        if(config.limitJsonOutput())
        {
            report = ReportLimiter.limitAllListsToMaxOne(report);
        }

        return report;
    }

    private Set<DoidNode> loadConfiguredPrimaryTumor(final OrangeConfig config)
    {
        Set<DoidNode> nodes = Sets.newHashSet();
        LOGGER.info("Determining configured primary tumor");
        for(String doid : config.primaryTumorDoids())
        {
            DoidNode node = resolveDoid(mDoidEntry.nodes(), doid);
            if(node != null)
            {
                LOGGER.info(" Adding DOID {} ({}) as configured primary tumor", doid, node.doidTerm());
                nodes.add(node);
            }
            else
            {
                LOGGER.warn("Could not resolve doid '{}'", doid);
            }
        }
        return nodes;
    }

    @Nullable
    private static DoidNode resolveDoid(final List<DoidNode> nodes, final String doid)
    {
        for(DoidNode node : nodes)
        {
            if(node.doid().equals(doid))
            {
                return node;
            }
        }
        return null;
    }

    @Nullable
    private static String determinePipelineVersion(final OrangeConfig config) throws IOException
    {
        String pipelineVersionFile = config.pipelineVersionFile();
        if(pipelineVersionFile == null)
        {
            LOGGER.warn("No pipeline version could be determined as pipeline version file was not passed");
            return null;
        }

        String pipelineVersion = PipelineVersionFile.majorDotMinorVersion(pipelineVersionFile);
        if(pipelineVersion != null)
        {
            LOGGER.info("Determined pipeline version to be 'v{}'", pipelineVersion);
        }
        else
        {
            LOGGER.warn("No pipeline version could be determined as version could not be resolved from {}", pipelineVersionFile);
        }
        return pipelineVersion;
    }

    @Nullable
    private static OrangeSample loadSampleData(final OrangeConfig config, boolean loadTumorSample) throws IOException
    {
        if(loadTumorSample)
        {
            LOGGER.info("Loading tumor sample data");
        }
        else
        {
            if(config.wgsRefConfig() != null && config.wgsRefConfig().refSampleWGSMetricsFile() != null
                    && config.wgsRefConfig().refSampleFlagstatFile() != null)
            {
                LOGGER.info("Loading reference sample data");
            }
            else
            {
                LOGGER.info("Skipping loading of reference sample data as no flagstat or WGS metrics has been provided");
                return null;
            }
        }

        String metricsFile = loadTumorSample ? config.tumorSampleWGSMetricsFile() : config.wgsRefConfig().refSampleWGSMetricsFile();
        WGSMetrics metrics = OrangeConversion.convert(BamMetricSummary.read(metricsFile));
        LOGGER.info(" Loaded WGS metrics from {}", metricsFile);

        String flagstatFile = loadTumorSample ? config.tumorSampleFlagstatFile() : config.wgsRefConfig().refSampleFlagstatFile();
        Flagstat flagstat = OrangeConversion.convert(BamFlagStats.read(flagstatFile));
        LOGGER.info(" Loaded flagstat from {}", flagstatFile);

        return ImmutableOrangeSample.builder().metrics(metrics).flagstat(flagstat).build();
    }

    @NotNull
    private static EnsemblDataCache loadEnsemblDataCache(final OrangeConfig config)
    {
        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(config.ensemblDataDirectory(),
                RefGenomeVersion.from(config.refGenomeVersion().name()));
        ensemblDataCache.setRequireNonEnsemblTranscripts();
        ensemblDataCache.load(false);
        return ensemblDataCache;
    }

    @Nullable
    private static Map<String, Double> loadGermlineMVLHPerGene(final OrangeConfig config, final Map<String,DriverGene> driverGenes)
            throws IOException
    {
        OrangeWGSRefConfig orangeWGSRefConfig = config.wgsRefConfig();
        String germlineGeneCoverageTsv = orangeWGSRefConfig != null ? orangeWGSRefConfig.germlineGeneCoverageTsv() : null;
        if(germlineGeneCoverageTsv == null)
        {
            LOGGER.info("Skipping loading of germline MVLH as no germline gene coverage has been provided");
            return null;
        }

        Map<String, Double> mvlhPerGene = GermlineMVLHFactory.loadGermlineMVLHPerGene(germlineGeneCoverageTsv, driverGenes);
        LOGGER.info("Loaded MVLH data for {} genes", mvlhPerGene.keySet().size());

        return mvlhPerGene;
    }

    private PurpleData loadPurpleData(final OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading PURPLE data from {}", config.purpleDataDirectory());

        String referenceSample = config.wgsRefConfig() != null ? config.wgsRefConfig().referenceSampleId() : null;

        PurpleData purple = PurpleDataLoader.load(config, mDriverGenes);

        DecimalFormat purityFormat = new DecimalFormat("#'%'");
        PurityContext purityContext = purple.purityContext();
        LOGGER.info("  QC status: {}", purityContext.qc().toString());
        LOGGER.info("  Tumor purity: {} ({}-{})",
                purityFormat.format(purityContext.bestFit().purity() * 100),
                purityFormat.format(purityContext.score().minPurity() * 100),
                purityFormat.format(purityContext.score().maxPurity() * 100));
        LOGGER.info("  Tumor ploidy: {} ({}-{})",
                purityContext.bestFit().ploidy(),
                purityContext.score().minPloidy(),
                purityContext.score().maxPloidy());
        LOGGER.info("  Fit method: {}", purityContext.method());
        LOGGER.info("  Whole genome duplication: {}", purityContext.wholeGenomeDuplication() ? "yes" : "no");
        LOGGER.info("  Microsatellite status: {}", purityContext.microsatelliteStatus().display());
        LOGGER.info("  Tumor mutational load status: {}", purityContext.tumorMutationalLoadStatus().display());

        LOGGER.info(" Loaded {} somatic driver catalog entries", purple.somaticDrivers().size());

        LOGGER.info(" Loaded {} somatic variants", purple.somaticVariants().size());

        LOGGER.info(" Loaded {} somatic copy numbers entries", purple.somaticCopyNumbers().size());
        LOGGER.info(" Loaded {} somatic gene copy numbers entries", purple.somaticGeneCopyNumbers().size());

        if(referenceSample != null)
        {
            LOGGER.info(" Loaded {} germline driver catalog entries", purple.germlineDrivers().size());
            LOGGER.info(" Loaded {} reportable germline variants", purple.germlineVariants().size());

            LOGGER.info(" Loaded {} reportable germline copy-number events", purple.germlineDeletions().size());
        }
        else
        {
            LOGGER.debug(" Skipped loading germline variants and deletions since no reference sample configured");
        }

        return purple;
    }

    private static LinxData loadLinxData(final OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading LINX somatic data from {}", config.linxSomaticDataDirectory());

        String linxGermlineDataDirectory = config.wgsRefConfig() != null ? config.wgsRefConfig().linxGermlineDataDirectory() : null;

        LinxData linx = LinxDataLoader.load(config);

        LOGGER.info(" Loaded {} somatic structural variants", linx.somaticSvAnnotations().size());
        LOGGER.info(" Loaded {} somatic structural drivers", linx.somaticDrivers().size());
        LOGGER.info(" Loaded {} somatic fusions", linx.fusions().size());

        LOGGER.info(" Loaded {} somatic breakends", linx.somaticBreakends().size());
        LOGGER.info(" Loaded {} somatic reportable homozygous disruptions", linx.somaticHomozygousDisruptions().size());

        if(linxGermlineDataDirectory != null)
        {
            LOGGER.info("Loading LINX germline data from {}", linxGermlineDataDirectory);
            LOGGER.info(" Loaded {} germline structural variants", linx.germlineSvAnnotations().size());
            LOGGER.info(" Loaded {} germline breakends (of which {} are reportable)",
                    linx.germlineBreakends().size(),
                    linx.germlineBreakends().size());
            LOGGER.info(" Loaded {} germline disruptions", linx.germlineDisruptions().size());
            LOGGER.info(" Loaded {} germline reportable homozygous disruptions", linx.germlineHomozygousDisruptions().size());
        }
        else
        {
            LOGGER.info(" Skipped loading LINX germline data as no linx germline data directory has been provided");
        }

        return linx;
    }

    @Nullable
    private IsofoxData loadIsofoxData(final OrangeConfig config) throws IOException
    {
        OrangeRnaConfig rna = config.rnaConfig();
        if(rna == null)
        {
            LOGGER.info("Skipping ISOFOX data loading as RNA is not configured");
            return null;
        }

        String orangeCancerType = mCohortMapper.cancerTypeForSample(createSample(config));
        if(orangeCancerType == null)
        {
            LOGGER.warn("Could not resolve ORANGE cancer type for {}", config.tumorSampleId());
            return null;
        }

        String isofoxCancerType;
        // TODO (KD): Replace with unified cohort mapping code (see also ACTIN-1010)
        if(orangeCancerType.equals("Ovary/Fallopian tube"))
        {
            isofoxCancerType = "Ovary";
            LOGGER.debug("Converted orange cancer type '{}' to isofox cancer type '{}'", orangeCancerType, isofoxCancerType);
        }
        else if(orangeCancerType.equals("Unknown"))
        {
            isofoxCancerType = null;
            LOGGER.debug("Converted orange cancer type '{}' to null isofox cancer type", orangeCancerType);
        }
        else
        {
            isofoxCancerType = orangeCancerType;
        }

        return IsofoxDataLoader.load(isofoxCancerType,
                rna.isofoxSummaryCsv(),
                rna.isofoxGeneDataCsv(),
                rna.isofoxFusionCsv(),
                rna.isofoxAltSpliceJunctionCsv());
    }

    @Nullable
    private static LilacSummaryData loadLilacData(final OrangeConfig config) throws IOException
    {
        if(config.lilacResultTsv() == null || config.lilacQcTsv() == null)
        {
            LOGGER.info("Skipping loading LILAC results since LILAC input dir or tsvs were not provided");
            return null;
        }

        return LilacSummaryData.load(config.lilacQcTsv(), config.lilacResultTsv());
    }

    @Nullable
    private static VirusInterpreterData loadVirusInterpreterData(final OrangeConfig config) throws IOException
    {
        if(config.wgsRefConfig() == null)
        {
            return null;
        }

        String annotatedVirusTsv = config.wgsRefConfig().annotatedVirusTsv();
        if(annotatedVirusTsv == null)
        {
            LOGGER.debug("Skipping loading of annotated viruses as no input has been provided");
            return null;
        }

        return VirusInterpreterDataLoader.load(annotatedVirusTsv);
    }

    @Nullable
    private static ChordData loadChordAnalysis(final OrangeConfig config) throws IOException
    {
        if(config.wgsRefConfig() == null)
        {
            return null;
        }

        String chordPredictionTxt = config.wgsRefConfig().chordPredictionTxt();
        if(chordPredictionTxt == null)
        {
            LOGGER.debug("Skipping CHORD loading as no input has been provided");
            return null;
        }

        LOGGER.info("Loading CHORD data from {}", new File(chordPredictionTxt).getParent());
        ChordData chordData = ChordDataFile.read(chordPredictionTxt);
        LOGGER.info(" HR Status: {} with type '{}'", chordData.hrStatus().display(), chordData.hrdType());
        return chordData;
    }

    @Nullable
    private static CuppaData loadCuppaData(final OrangeConfig config) throws Exception
    {
        OrangeWGSRefConfig orangeWGSRefConfig = config.wgsRefConfig();
        if(orangeWGSRefConfig == null || orangeWGSRefConfig.cuppaVisDataTsv() == null)
        {
            return null;
        }

        String cuppaVisDataTsv = orangeWGSRefConfig.cuppaVisDataTsv();
        LOGGER.info("Loading CUPPA predictions from {}", new File(cuppaVisDataTsv).getParent());
        CuppaData cuppaData = CuppaDataFactory.create(cuppaVisDataTsv);
        LOGGER.info(" Loaded {} CUPPA predictions from {}", cuppaData.predictions().size(), cuppaVisDataTsv);

        return cuppaData;
    }

    @Nullable
    private static List<PeachGenotype> loadPeachData(final OrangeConfig config) throws IOException
    {
        String peachGenotypeTsv = config.wgsRefConfig() != null ? config.wgsRefConfig().peachGenotypeTsv() : null;

        if(peachGenotypeTsv == null)
        {
            LOGGER.info("Skipping PEACH loading since no peach genotype tsv has been provided");
            return null;
        }

        LOGGER.info("Loading PEACH from {}", new File(peachGenotypeTsv).getParent());
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(peachGenotypeTsv);
        LOGGER.info(" Loaded {} PEACH genotypes from {}", peachGenotypes.size(), config.wgsRefConfig().peachGenotypeTsv());
        List<PeachGenotype> filterUGT1A1FromPeachGenotypes = peachGenotypes.stream().filter(genotype -> !genotype.gene().equals("UGT1A1")).toList();
        return filterUGT1A1FromPeachGenotypes;
    }

    @Nullable
    private static List<SignatureAllocation> loadSigAllocations(final OrangeConfig config) throws IOException
    {
        if(config.wgsRefConfig() == null)
        {
            return null;
        }

        String sigsAllocationTsv = config.wgsRefConfig().sigsAllocationTsv();

        if(sigsAllocationTsv == null)
        {
            LOGGER.info("Skipping signature loading since no sigs allocation tsv has been provided");
            return null;
        }

        LOGGER.info("Loading Sigs from {}", new File(sigsAllocationTsv).getParent());
        List<SignatureAllocation> sigsAllocations = SignatureAllocationFile.read(sigsAllocationTsv);
        LOGGER.info(" Loaded {} signature allocations from {}", sigsAllocations.size(), sigsAllocationTsv);

        return sigsAllocations;
    }

    private OrangePlots buildPlots(final OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading plots");

        mPlotManager.createPlotDirectory();

        String linxPlotDir = config.linxPlotDirectory();
        List<String> linxDriverPlots = Lists.newArrayList();

        if(linxPlotDir != null)
        {
            for(String file : new File(linxPlotDir).list())
            {
                linxDriverPlots.add(mPlotManager.processPlotFile(linxPlotDir + File.separator + file));
            }

            LOGGER.info(" Loaded {} linx plots from {}", linxDriverPlots.size(), linxPlotDir);
        }

        String referenceBqrPlot = mPlotManager.processPlotFile(
                config.wgsRefConfig() != null ? config.wgsRefConfig().refSampleBqrPlot() : null);

        String tumorBqrPlot = mPlotManager.processPlotFile(config.tumorSampleBqrPlot());

        String purplePlotBasePath = config.purplePlotDirectory() + File.separator + config.tumorSampleId();
        String purpleInputPlot = mPlotManager.processPlotFile(purplePlotBasePath + ".input.png");
        String purpleFinalCircosPlot = mPlotManager.processPlotFile(purplePlotBasePath + ".circos.png");
        String purpleClonalityPlot = mPlotManager.processPlotFile(purplePlotBasePath + ".somatic.clonality.png");
        String purpleCopyNumberPlot = mPlotManager.processPlotFile(purplePlotBasePath + ".copynumber.png");
        String purpleVariantCopyNumberPlot = mPlotManager.processPlotFile(purplePlotBasePath + ".somatic.png");
        String purplePurityRangePlot = mPlotManager.processPlotFile(purplePlotBasePath + ".purity.range.png");
        String purpleKataegisPlot = mPlotManager.processPlotFile(purplePlotBasePath + ".somatic.rainfall.png");

        List<String> purplePlots = Arrays.asList(purpleInputPlot, purpleFinalCircosPlot, purpleClonalityPlot, purpleCopyNumberPlot,
                purpleVariantCopyNumberPlot, purplePurityRangePlot, purpleKataegisPlot);

        if(purplePlots.stream().anyMatch(Objects::isNull))
        {
            LOGGER.warn("Skipping making ORANGE report: missing one or more PURPLE plot paths, likely because the input sample(s) has no or extremely sparse data");
            System.exit(0);
        }

        String cuppaSummaryPlot = mPlotManager.processPlotFile(
                (config.wgsRefConfig() != null) ? config.wgsRefConfig().cuppaSummaryPlot() : null
        );

        return ImmutableOrangePlots.builder()
                .sageReferenceBQRPlot(referenceBqrPlot)
                .sageTumorBQRPlot(tumorBqrPlot)
                .purpleInputPlot(purpleInputPlot)
                .purpleFinalCircosPlot(purpleFinalCircosPlot)
                .purpleClonalityPlot(purpleClonalityPlot)
                .purpleCopyNumberPlot(purpleCopyNumberPlot)
                .purpleVariantCopyNumberPlot(purpleVariantCopyNumberPlot)
                .purplePurityRangePlot(purplePurityRangePlot)
                .purpleKataegisPlot(purpleKataegisPlot)
                .linxDriverPlots(linxDriverPlots)
                .cuppaSummaryPlot(cuppaSummaryPlot)
                .build();
    }

    private static void verifyPlots(final OrangePlots orangePlots, final LinxData linxData)
    {
        Set<Integer> linxVisualizedClusters = LinxReportableClusters.findVisualizedClusters(linxData);

        if(linxVisualizedClusters.size() != orangePlots.linxDriverPlots().size())
        {
            LOGGER.warn("Expected {} linx plots, but found {}", linxVisualizedClusters.size(), orangePlots.linxDriverPlots().size());
        }
    }

    private static Sample createSample(final OrangeConfig config)
    {
        return ImmutableSample.builder()
                .sampleId(config.tumorSampleId())
                .doids(config.primaryTumorDoids())
                .build();
    }

    @VisibleForTesting
    public void setSuppressGeneWarnings()
    {
        mSuppressGeneWarnings = true;
    }
}
