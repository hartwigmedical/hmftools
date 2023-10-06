package com.hartwig.hmftools.orange.algo;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.isofox.IsofoxData;
import com.hartwig.hmftools.common.isofox.IsofoxDataLoader;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDataLoader;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.sage.GeneDepthFile;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.datamodel.cohort.Evaluation;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.flagstat.Flagstat;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.metrics.WGSMetrics;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangePlots;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeSample;
import com.hartwig.hmftools.datamodel.orange.OrangePlots;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeSample;
import com.hartwig.hmftools.datamodel.orange.PercentileType;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.wildtype.WildTypeGene;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeRnaConfig;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaDataFactory;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpreter;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpreter;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.plot.DummyPlotManager;
import com.hartwig.hmftools.orange.algo.plot.FileBasedPlotManager;
import com.hartwig.hmftools.orange.algo.plot.PlotManager;
import com.hartwig.hmftools.orange.algo.purple.GermlineGainLossFactory;
import com.hartwig.hmftools.orange.algo.purple.PurpleData;
import com.hartwig.hmftools.orange.algo.purple.PurpleDataLoader;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpreter;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariantFactory;
import com.hartwig.hmftools.orange.algo.util.GermlineConversion;
import com.hartwig.hmftools.orange.algo.util.ReportLimiter;
import com.hartwig.hmftools.orange.algo.wildtype.WildTypeAlgo;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableObservation;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableSample;
import com.hartwig.hmftools.orange.cohort.datamodel.Observation;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapper;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapping;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMappingFile;
import com.hartwig.hmftools.orange.cohort.mapping.DoidCohortMapper;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentiles;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentilesFile;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentilesModel;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.OrangeConversion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class OrangeAlgo
{
    @NotNull
    private final DoidEntry doidEntry;
    @NotNull
    private final CohortMapper cohortMapper;
    @NotNull
    private final CohortPercentilesModel percentilesModel;
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;
    @NotNull
    private final EnsemblDataCache ensemblDataCache;
    @NotNull
    private final PlotManager plotManager;

    private boolean suppressGeneWarnings;

    @NotNull
    public static OrangeAlgo fromConfig(@NotNull OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading DOID database from {}", config.doidJsonFile());
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJsonFile());
        DoidParents doidParentModel = DoidParents.fromEdges(doidEntry.edges());

        LOGGER.info("Reading cohort mappings from {}", config.cohortMappingTsv());
        List<CohortMapping> mappings = CohortMappingFile.read(config.cohortMappingTsv());
        LOGGER.info(" Reading {} cohort mappings", mappings.size());
        CohortMapper mapper = new DoidCohortMapper(doidParentModel, mappings);

        LOGGER.info("Reading percentiles from {}", config.cohortPercentilesTsv());
        Multimap<PercentileType, CohortPercentiles> percentilesMap = CohortPercentilesFile.read(config.cohortPercentilesTsv());
        LOGGER.info(" Read {} percentiles", percentilesMap.values().size());
        CohortPercentilesModel percentilesModel = new CohortPercentilesModel(mapper, percentilesMap);

        LOGGER.info("Reading driver genes from {}", config.driverGenePanelTsv());
        List<DriverGene> driverGenes = DriverGeneFile.read(config.driverGenePanelTsv());
        LOGGER.info(" Read {} driver genes", driverGenes.size());

        LOGGER.info("Reading known fusions from {}", config.knownFusionFile());
        KnownFusionCache knownFusionCache = new KnownFusionCache();
        if(!knownFusionCache.loadFile(config.knownFusionFile()))
        {
            throw new IOException("Could not load known fusions from " + config.knownFusionFile());
        }
        LOGGER.info(" Read {} known fusion entries", knownFusionCache.getData().size());

        LOGGER.info("Reading ensembl data cache from {}", config.ensemblDataDirectory());
        EnsemblDataCache ensemblDataCache = loadEnsemblDataCache(config);
        LOGGER.info(" Read ensembl data dir");

        String outputDir = config.outputDir();
        PlotManager plotManager = !outputDir.isEmpty() ? new FileBasedPlotManager(outputDir) : new DummyPlotManager();

        return new OrangeAlgo(doidEntry, mapper, percentilesModel, driverGenes, knownFusionCache, ensemblDataCache, plotManager);
    }

    private OrangeAlgo(@NotNull final DoidEntry doidEntry, @NotNull final CohortMapper cohortMapper,
            @NotNull final CohortPercentilesModel percentilesModel, @NotNull final List<DriverGene> driverGenes,
            @NotNull final KnownFusionCache knownFusionCache, @NotNull final EnsemblDataCache ensemblDataCache,
            @NotNull final PlotManager plotManager)
    {
        this.doidEntry = doidEntry;
        this.cohortMapper = cohortMapper;
        this.percentilesModel = percentilesModel;
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
        this.ensemblDataCache = ensemblDataCache;
        this.plotManager = plotManager;
        suppressGeneWarnings = false;
    }

    @NotNull
    public OrangeRecord run(@NotNull OrangeConfig config) throws IOException
    {
        Set<DoidNode> configuredPrimaryTumor = loadConfiguredPrimaryTumor(config);
        String platinumVersion = determinePlatinumVersion(config);
        OrangeSample refSample = loadSampleData(config, false);
        OrangeSample tumorSample = loadSampleData(config, true);

        PurpleData purpleData = loadPurpleData(config);
        LinxData linxData = loadLinxData(config);
        Map<String, Double> mvlhPerGene = loadGermlineMVLHPerGene(config);
        ChordData chord = loadChordAnalysis(config);
        LilacSummaryData lilac = loadLilacData(config);
        VirusInterpreterData virusInterpreter = loadVirusInterpreterData(config);
        CuppaData cuppa = loadCuppaData(config);
        List<PeachGenotype> peach = loadPeachData(config);
        List<SignatureAllocation> sigAllocations = loadSigAllocations(config);
        IsofoxData isofoxData = loadIsofoxData(config);

        ExperimentType experimentType = purpleData.purityContext().targeted() ? ExperimentType.TARGETED : ExperimentType.WHOLE_GENOME;
        LOGGER.info("Determined experiment type to be '{}'", experimentType);

        LinxInterpreter linxInterpreter = new LinxInterpreter(driverGenes, knownFusionCache);
        LinxRecord linx = linxInterpreter.interpret(linxData);

        PaveAlgo pave = new PaveAlgo(ensemblDataCache, !suppressGeneWarnings);

        PurpleVariantFactory purpleVariantFactory = new PurpleVariantFactory(pave);
        GermlineGainLossFactory germlineGainLossFactory = new GermlineGainLossFactory(ensemblDataCache);
        PurpleInterpreter purpleInterpreter =
                new PurpleInterpreter(purpleVariantFactory, germlineGainLossFactory, driverGenes, linx, chord);
        PurpleRecord purple = purpleInterpreter.interpret(purpleData);

        IsofoxRecord isofox = null;
        if(isofoxData != null)
        {
            IsofoxInterpreter isofoxInterpreter = new IsofoxInterpreter(driverGenes, knownFusionCache, linx);
            isofox = isofoxInterpreter.interpret(isofoxData);
        }

        List<WildTypeGene> wildTypeGenes = Lists.newArrayList();
        if(WildTypeAlgo.wildTypeCallingAllowed(purple.fit().qc().status()))
        {
            wildTypeGenes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                    purple.reportableSomaticVariants(),
                    purple.reportableGermlineVariants(),
                    purple.reportableSomaticGainsLosses(),
                    linx.reportableSomaticFusions(),
                    linx.somaticHomozygousDisruptions(),
                    linx.reportableSomaticBreakends());
            LOGGER.info("Identified {} of {} driver genes to be wild-type", wildTypeGenes.size(), driverGenes.size());
        }
        else
        {
            LOGGER.info("Wild-type calling skipped due to insufficient tumor sample quality");
        }

        OrangeRecord report = ImmutableOrangeRecord.builder()
                .sampleId(config.tumorSampleId())
                .samplingDate(config.experimentDate())
                .experimentType(experimentType)
                .configuredPrimaryTumor(ConversionUtil.mapToIterable(configuredPrimaryTumor, OrangeConversion::convert))
                .refGenomeVersion(config.refGenomeVersion())
                .platinumVersion(platinumVersion)
                .refSample(refSample)
                .tumorSample(tumorSample)
                .germlineMVLHPerGene(mvlhPerGene)
                .purple(purple)
                .linx(linx)
                .wildTypeGenes(wildTypeGenes)
                .isofox(isofox)
                .lilac(OrangeConversion.convert(lilac))
                .virusInterpreter(virusInterpreter != null ? OrangeConversion.convert(virusInterpreter) : null)
                .chord(chord != null ? OrangeConversion.convert(chord) : null)
                .cuppa(cuppa)
                .peach(ConversionUtil.mapToIterable(peach, OrangeConversion::convert))
                .sigAllocations(ConversionUtil.mapToIterable(sigAllocations, OrangeConversion::convert))
                .cohortEvaluations(evaluateCohortPercentiles(config, purple))
                .plots(buildPlots(config))
                .build();

        if(config.limitJsonOutput())
        {
            report = ReportLimiter.limitAllListsToMaxOne(report);
        }

        if(config.convertGermlineToSomatic())
        {
            report = GermlineConversion.convertGermlineToSomatic(report);
        }

        return report;
    }

    @NotNull
    private Set<DoidNode> loadConfiguredPrimaryTumor(@NotNull OrangeConfig config)
    {
        Set<DoidNode> nodes = Sets.newHashSet();
        LOGGER.info("Determining configured primary tumor");
        for(String doid : config.primaryTumorDoids())
        {
            DoidNode node = resolveDoid(doidEntry.nodes(), doid);
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
    private static DoidNode resolveDoid(@NotNull List<DoidNode> nodes, @NotNull String doid)
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
    private static String determinePlatinumVersion(@NotNull OrangeConfig config) throws IOException
    {
        String pipelineVersionFile = config.pipelineVersionFile();
        if(pipelineVersionFile == null)
        {
            LOGGER.warn("No platinum version could be determined as pipeline version file was not passed");
            return null;
        }

        String platinumVersion = PipelineVersionFile.majorDotMinorVersion(pipelineVersionFile);
        if(platinumVersion != null)
        {
            LOGGER.info("Determined platinum version to be 'v{}'", platinumVersion);
        }
        else
        {
            LOGGER.warn("No platinum version could be determined as version could not be resolved from {}", pipelineVersionFile);
        }
        return platinumVersion;
    }

    @Nullable
    private static OrangeSample loadSampleData(@NotNull OrangeConfig config, boolean loadTumorSample) throws IOException
    {
        if(loadTumorSample)
        {
            LOGGER.info("Loading tumor sample data");
        }
        else
        {
            if(config.refSampleWGSMetricsFile() != null && config.refSampleFlagstatFile() != null)
            {
                LOGGER.info("Loading reference sample data");
            }
            else
            {
                LOGGER.info("Skipping loading of reference sample data as no flagstat or WGS metrics has been provided");
                return null;
            }
        }

        String metricsFile = loadTumorSample ? config.tumorSampleWGSMetricsFile() : config.refSampleWGSMetricsFile();
        WGSMetrics metrics = OrangeConversion.convert(WGSMetricsFile.read(metricsFile));
        LOGGER.info(" Loaded WGS metrics from {}", metricsFile);

        String flagstatFile = loadTumorSample ? config.tumorSampleFlagstatFile() : config.refSampleFlagstatFile();
        Flagstat flagstat = OrangeConversion.convert(FlagstatFile.read(flagstatFile));
        LOGGER.info(" Loaded flagstat from {}", flagstatFile);

        return ImmutableOrangeSample.builder().metrics(metrics).flagstat(flagstat).build();
    }

    @NotNull
    private static EnsemblDataCache loadEnsemblDataCache(@NotNull OrangeConfig config)
    {
        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(config.ensemblDataDirectory(),
                RefGenomeVersion.from(config.refGenomeVersion().name()));
        ensemblDataCache.setRequireNonEnsemblTranscripts();
        ensemblDataCache.load(false);
        return ensemblDataCache;
    }

    @Nullable
    private static Map<String, Double> loadGermlineMVLHPerGene(@NotNull OrangeConfig config) throws IOException
    {
        String sageGermlineGeneCoverageTsv = config.sageGermlineGeneCoverageTsv();
        if(sageGermlineGeneCoverageTsv == null)
        {
            LOGGER.info("Skipping loading of germline MVLH as no germline gene coverage has been provided");
            return null;
        }

        Map<String, Double> mvlhPerGene = Maps.newTreeMap();

        List<String> lines = Files.readAllLines(new File(sageGermlineGeneCoverageTsv).toPath());
        String header = lines.get(0);

        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        int geneIndex = fieldsIndexMap.get(GeneDepthFile.COL_GENE);
        int mvlhIndex = fieldsIndexMap.get(GeneDepthFile.COL_MV_LIKELIHOOD);

        for(String line : lines.subList(1, lines.size()))
        {
            String[] values = line.split(TSV_DELIM);
            String gene = values[geneIndex];
            String mvlhString = values[mvlhIndex].substring(0, values[mvlhIndex].length() - 1);
            double missedVariantLikelihood = Double.parseDouble(mvlhString) / 100D;
            mvlhPerGene.put(gene, missedVariantLikelihood);
        }

        LOGGER.info("Loaded MVLH data for {} genes", mvlhPerGene.keySet().size());
        return mvlhPerGene;
    }

    @NotNull
    private static PurpleData loadPurpleData(@NotNull OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading PURPLE data from {}", config.purpleDataDirectory());

        String referenceSample = config.referenceSampleId();
        PurpleData purple = PurpleDataLoader.load(config.tumorSampleId(),
                referenceSample,
                config.rnaConfig() != null ? config.rnaConfig().rnaSampleId() : null,
                config.purpleDataDirectory());

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

        LOGGER.info(" Loaded {} somatic variants (of which {} are reportable)",
                purple.allSomaticVariants().size(),
                purple.reportableSomaticVariants().size());
        LOGGER.info(" Loaded {} somatic copy numbers entries", purple.allSomaticCopyNumbers().size());
        LOGGER.info(" Loaded {} somatic gene copy numbers entries", purple.allSomaticGeneCopyNumbers().size());
        LOGGER.info(" Loaded {} somatic structural variants", purple.allSomaticStructuralVariants().size());

        if(referenceSample != null)
        {
            LOGGER.info(" Loaded {} germline driver catalog entries", purple.germlineDrivers().size());
            LOGGER.info(" Loaded {} germline variants (of which {} are reportable)",
                    purple.allGermlineVariants().size(),
                    purple.reportableGermlineVariants().size());

            LOGGER.info(" Loaded {} germline deletions (of which {} are reportable)",
                    purple.allGermlineDeletions().size(),
                    purple.reportableGermlineDeletions().size());

            LOGGER.info(" Loaded {} germline structural variants", purple.allGermlineStructuralVariants().size());
        }
        else
        {
            LOGGER.debug(" Skipped loading germline variants and deletions since no reference sample configured");
        }

        return purple;
    }

    @NotNull
    private static LinxData loadLinxData(@NotNull OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading LINX somatic data from {}", config.linxSomaticDataDirectory());

        String linxGermlineDataDirectory = config.linxGermlineDataDirectory();
        LinxData linx = LinxDataLoader.load(config.tumorSampleId(), config.linxSomaticDataDirectory(), linxGermlineDataDirectory);

        LOGGER.info(" Loaded {} somatic structural variants", linx.allSomaticStructuralVariants().size());
        LOGGER.info(" Loaded {} somatic structural drivers", linx.somaticDrivers().size());
        LOGGER.info(" Loaded {} somatic fusions (of which {} are reportable)",
                linx.allSomaticFusions().size(),
                linx.reportableSomaticFusions().size());
        LOGGER.info(" Loaded {} somatic breakends (of which {} are reportable)",
                linx.allSomaticBreakends().size(),
                linx.reportableSomaticBreakends().size());
        LOGGER.info(" Loaded {} somatic reportable homozygous disruptions", linx.somaticHomozygousDisruptions().size());

        if(!config.tumorOnlyMode() && linxGermlineDataDirectory != null)
        {
            LOGGER.info("Loading LINX germline data from {}", linxGermlineDataDirectory);
            LOGGER.info(" Loaded {} germline structural variants", linx.allGermlineStructuralVariants().size());
            LOGGER.info(" Loaded {} germline breakends (of which {} are reportable)",
                    linx.allGermlineBreakends().size(),
                    linx.reportableGermlineBreakends().size());
            LOGGER.info(" Loaded {} germline disruptions (of which {} are reportable)",
                    linx.allGermlineDisruptions().size(),
                    linx.reportableGermlineDisruptions().size());
            LOGGER.info(" Loaded {} germline reportable homozygous disruptions", linx.germlineHomozygousDisruptions().size());
        }
        else
        {
            LOGGER.info(" Skipped loading LINX germline data as no linx germline data directory has been provided");
        }

        return linx;
    }

    @Nullable
    private IsofoxData loadIsofoxData(@NotNull OrangeConfig config) throws IOException
    {
        OrangeRnaConfig rna = config.rnaConfig();
        if(rna == null)
        {
            LOGGER.info("Skipping ISOFOX data loading as RNA is not configured");
            return null;
        }

        String isofoxCancerType = cohortMapper.cancerTypeForSample(createSample(config));
        if(isofoxCancerType == null)
        {
            LOGGER.warn("Could not resolve isofox cancer type for {}" + config.tumorSampleId());
            return null;
        }

        return IsofoxDataLoader.load(isofoxCancerType,
                rna.isofoxGeneDistributionCsv(),
                rna.isofoxAltSjCohortCsv(),
                rna.isofoxSummaryCsv(),
                rna.isofoxGeneDataCsv(),
                rna.isofoxFusionCsv(),
                rna.isofoxAltSpliceJunctionCsv());
    }

    @NotNull
    private static LilacSummaryData loadLilacData(@NotNull OrangeConfig config) throws IOException
    {
        return LilacSummaryData.load(config.lilacQcCsv(), config.lilacResultCsv());
    }

    @Nullable
    private static VirusInterpreterData loadVirusInterpreterData(@NotNull OrangeConfig config) throws IOException
    {
        if(config.tumorOnlyMode())
        {
            return null;
        }

        String annotatedVirusTsv = config.annotatedVirusTsv();
        if(annotatedVirusTsv == null)
        {
            LOGGER.debug("Skipping loading of annotated  viruses as no input has been provided");
            return null;
        }

        return VirusInterpreterDataLoader.load(annotatedVirusTsv);
    }

    @Nullable
    private static ChordData loadChordAnalysis(@NotNull OrangeConfig config) throws IOException
    {
        if(config.tumorOnlyMode())
        {
            return null;
        }

        String chordPredictionTxt = config.chordPredictionTxt();
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
    private static CuppaData loadCuppaData(@NotNull OrangeConfig config) throws IOException
    {
        if(config.tumorOnlyMode())
        {
            return null;
        }

        String cuppaResultTsv = config.cuppaResultCsv();
        if(cuppaResultTsv == null)
        {
            LOGGER.debug("Skipping CUPPA loading as no input has been provided");
            return null;
        }

        LOGGER.info("Loading CUPPA from {}", new File(cuppaResultTsv).getParent());
        List<CuppaDataFile> cuppaEntries = CuppaDataFile.read(cuppaResultTsv);
        LOGGER.info(" Loaded {} entries from {}", cuppaEntries.size(), cuppaResultTsv);

        CuppaData cuppaData = CuppaDataFactory.create(cuppaEntries);
        CuppaPrediction best = cuppaData.predictions().get(0);
        LOGGER.info(" Predicted cancer type '{}' with likelihood {}", best.cancerType(), best.likelihood());

        return cuppaData;
    }

    @Nullable
    private static List<PeachGenotype> loadPeachData(@NotNull OrangeConfig config) throws IOException
    {
        String peachGenotypeTsv = config.peachGenotypeTsv();

        if(peachGenotypeTsv == null)
        {
            LOGGER.info("Skipping PEACH loading since no peach genotype tsv has been provided");
            return null;
        }

        LOGGER.info("Loading PEACH from {}", new File(peachGenotypeTsv).getParent());
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(peachGenotypeTsv);
        LOGGER.info(" Loaded {} PEACH genotypes from {}", peachGenotypes.size(), config.peachGenotypeTsv());

        return peachGenotypes;
    }

    @Nullable
    private static List<SignatureAllocation> loadSigAllocations(@NotNull OrangeConfig config) throws IOException
    {
        if(config.tumorOnlyMode())
        {
            return null;
        }

        String sigsAllocationTsv = config.sigsAllocationTsv();

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

    @NotNull
    private Map<PercentileType, Evaluation> evaluateCohortPercentiles(@NotNull OrangeConfig config, @NotNull PurpleRecord purple)
    {
        PercentileType type = PercentileType.SV_TMB;

        Observation svTmbObservation = ImmutableObservation.builder()
                .sample(createSample(config))
                .type(type)
                .value(purple.characteristics().svTumorMutationalBurden())
                .build();

        LOGGER.info("Determining SV TMB percentile for value {}", svTmbObservation.value());
        Map<PercentileType, Evaluation> evaluations = Maps.newHashMap();
        Evaluation evaluation = percentilesModel.percentile(svTmbObservation);
        if(evaluation != null)
        {
            String cancerType = evaluation.cancerType();
            Double cancerTypePercentile = evaluation.cancerTypePercentile();
            LOGGER.info(" Determined percentile '{}' for pan-cancer and '{}' for cancer type '{}'",
                    evaluation.panCancerPercentile(),
                    cancerTypePercentile != null ? cancerTypePercentile : "NA",
                    cancerType != null ? cancerType : "NA");
            evaluations.put(type, evaluation);
        }
        else
        {
            LOGGER.warn("Could not evaluate SV TMB percentile for {}!", config.tumorSampleId());
        }

        return evaluations;
    }

    @NotNull
    private OrangePlots buildPlots(@NotNull OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading plots");

        plotManager.createPlotDirectory();

        String linxPlotDir = config.linxPlotDirectory();
        List<String> linxDriverPlots = Lists.newArrayList();
        if(new File(linxPlotDir).exists())
        {
            for(String file : new File(linxPlotDir).list())
            {
                linxDriverPlots.add(plotManager.processPlotFile(linxPlotDir + File.separator + file));
            }
            LOGGER.info(" Loaded {} linx plots from {}", linxDriverPlots.size(), linxPlotDir);
        }
        else
        {
            LOGGER.debug(" No linx plots have been loaded as plot directory {} does not exist", linxPlotDir);
        }

        String sageReferenceBQRPlot = plotManager.processPlotFile(config.sageSomaticRefSampleBQRPlot());
        String sageTumorBQRPlot = plotManager.processPlotFile(config.sageSomaticTumorSampleBQRPlot());

        String purplePlotBasePath = config.purplePlotDirectory() + File.separator + config.tumorSampleId();
        String purpleInputPlot = plotManager.processPlotFile(purplePlotBasePath + ".input.png");
        String purpleFinalCircosPlot = plotManager.processPlotFile(purplePlotBasePath + ".circos.png");
        String purpleClonalityPlot = plotManager.processPlotFile(purplePlotBasePath + ".somatic.clonality.png");
        String purpleCopyNumberPlot = plotManager.processPlotFile(purplePlotBasePath + ".copynumber.png");
        String purpleVariantCopyNumberPlot = plotManager.processPlotFile(purplePlotBasePath + ".somatic.png");
        String purplePurityRangePlot = plotManager.processPlotFile(purplePlotBasePath + ".purity.range.png");
        String purpleKataegisPlot = plotManager.processPlotFile(purplePlotBasePath + ".somatic.rainfall.png");

        String cuppaSummaryPlot = null;
        if(config.cuppaSummaryPlot() != null)
        {
            cuppaSummaryPlot = plotManager.processPlotFile(config.cuppaSummaryPlot());
        }

        String cuppaFeaturePlot = null;
        if(config.cuppaFeaturePlot() != null && new File(config.cuppaFeaturePlot()).exists())
        {
            cuppaFeaturePlot = plotManager.processPlotFile(config.cuppaFeaturePlot());
        }

        String cuppaChartPlot = null;
        if(config.cuppaChartPlot() != null)
        {
            cuppaChartPlot = plotManager.processPlotFile(config.cuppaChartPlot());
        }

        return ImmutableOrangePlots.builder()
                .sageReferenceBQRPlot(sageReferenceBQRPlot)
                .sageTumorBQRPlot(sageTumorBQRPlot)
                .purpleInputPlot(purpleInputPlot)
                .purpleFinalCircosPlot(purpleFinalCircosPlot)
                .purpleClonalityPlot(purpleClonalityPlot)
                .purpleCopyNumberPlot(purpleCopyNumberPlot)
                .purpleVariantCopyNumberPlot(purpleVariantCopyNumberPlot)
                .purplePurityRangePlot(purplePurityRangePlot)
                .purpleKataegisPlot(purpleKataegisPlot)
                .linxDriverPlots(linxDriverPlots)
                .cuppaSummaryPlot(cuppaSummaryPlot)
                .cuppaFeaturePlot(cuppaFeaturePlot)
                .cuppaChartPlot(cuppaChartPlot)
                .build();
    }

    @NotNull
    private static Sample createSample(@NotNull OrangeConfig config)
    {
        return ImmutableSample.builder().sampleId(config.tumorSampleId()).doids(config.primaryTumorDoids()).build();
    }

    @VisibleForTesting
    public void setSuppressGeneWarnings()
    {
        suppressGeneWarnings = true;
    }
}
