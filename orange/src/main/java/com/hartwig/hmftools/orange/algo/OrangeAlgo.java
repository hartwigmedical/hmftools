package com.hartwig.hmftools.orange.algo;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.isofox.IsofoxDataLoader;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDataLoader;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleDataLoader;
import com.hartwig.hmftools.common.sage.GeneDepthFile;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeRNAConfig;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaData;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaDataFactory;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaPrediction;
import com.hartwig.hmftools.orange.algo.interpretation.GermlineConversion;
import com.hartwig.hmftools.orange.algo.interpretation.ReportLimiter;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpretedData;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpreter;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpretedData;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpreter;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpreter;
import com.hartwig.hmftools.orange.algo.wildtype.WildTypeFactory;
import com.hartwig.hmftools.orange.algo.wildtype.WildTypeGene;
import com.hartwig.hmftools.orange.cohort.datamodel.Evaluation;
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
import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class OrangeAlgo {

    private static final Logger LOGGER = LogManager.getLogger(OrangeAlgo.class);

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
    public static OrangeAlgo fromConfig(@NotNull OrangeConfig config) throws IOException {
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
        if (!knownFusionCache.loadFile(config.knownFusionFile())) {
            throw new IOException("Could not load known fusions from " + config.knownFusionFile());
        }
        LOGGER.info(" Read {} known fusion entries", knownFusionCache.getData().size());

        return new OrangeAlgo(doidEntry, mapper, percentilesModel, driverGenes, knownFusionCache);
    }

    private OrangeAlgo(@NotNull final DoidEntry doidEntry, @NotNull final CohortMapper cohortMapper,
            @NotNull final CohortPercentilesModel percentilesModel, @NotNull final List<DriverGene> driverGenes,
            @NotNull final KnownFusionCache knownFusionCache) {
        this.doidEntry = doidEntry;
        this.cohortMapper = cohortMapper;
        this.percentilesModel = percentilesModel;
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
    }

    @NotNull
    public OrangeReport run(@NotNull OrangeConfig config) throws IOException {
        Set<DoidNode> configuredPrimaryTumor = loadConfiguredPrimaryTumor(config);
        String platinumVersion = determinePlatinumVersion(config);
        OrangeSample refSample = loadSampleData(config, false);
        OrangeSample tumorSample = loadSampleData(config, true);

        LinxInterpreter linxInterpreter = new LinxInterpreter(driverGenes, knownFusionCache);
        LinxInterpretedData linx = linxInterpreter.interpret(loadLinxData(config));

        ChordData chord = loadChordAnalysis(config);
        PurpleInterpreter purpleInterpreter = new PurpleInterpreter(driverGenes, chord);
        PurpleInterpretedData purple = purpleInterpreter.interpret(loadPurpleData(config));

        List<WildTypeGene> wildTypeGenes = WildTypeFactory.filterQCWildTypes(purple.fit().qc().status(),
                WildTypeFactory.determineWildTypeGenes(driverGenes,
                        purple.reportableSomaticVariants(),
                        purple.reportableGermlineVariants(),
                        purple.reportableSomaticGainsLosses(),
                        linx.reportableFusions(),
                        linx.homozygousDisruptions(),
                        linx.reportableBreakends()));
        LOGGER.info("Identified {} of {} driver genes to be wild-type", wildTypeGenes.size(), driverGenes.size());

        OrangeReport report = ImmutableOrangeReport.builder()
                .sampleId(config.tumorSampleId())
                .experimentDate(config.experimentDate())
                .configuredPrimaryTumor(configuredPrimaryTumor)
                .refGenomeVersion(config.refGenomeVersion())
                .platinumVersion(platinumVersion)
                .refSample(refSample)
                .tumorSample(tumorSample)
                .germlineMVLHPerGene(loadGermlineMVLHPerGene(config))
                .purple(purple)
                .linx(linx)
                .wildTypeGenes(wildTypeGenes)
                .isofox(loadIsofoxData(config, linx))
                .lilac(loadLilacData(config))
                .virusInterpreter(loadVirusInterpreterData(config))
                .chord(chord)
                .cuppa(loadCuppaData(config))
                .peach(loadPeachData(config))
                .cohortEvaluations(evaluateCohortPercentiles(config, purple))
                .plots(buildPlots(config))
                .build();

        if (config.limitJsonOutput()) {
            report = ReportLimiter.limitAllListsToMaxOne(report);
        }

        if (config.convertGermlineToSomatic()) {
            report = GermlineConversion.convertGermlineToSomatic(report);
        }

        return report;
    }

    @NotNull
    private Set<DoidNode> loadConfiguredPrimaryTumor(@NotNull OrangeConfig config) {
        Set<DoidNode> nodes = Sets.newHashSet();
        LOGGER.info("Determining configured primary tumor");
        for (String doid : config.primaryTumorDoids()) {
            DoidNode node = resolveDoid(doidEntry.nodes(), doid);
            if (node != null) {
                LOGGER.info(" Adding DOID {} ({}) as configured primary tumor", doid, node.doidTerm());
                nodes.add(node);
            } else {
                LOGGER.warn("Could not resolve doid '{}'", doid);
            }
        }
        return nodes;
    }

    @Nullable
    private static DoidNode resolveDoid(@NotNull List<DoidNode> nodes, @NotNull String doid) {
        for (DoidNode node : nodes) {
            if (node.doid().equals(doid)) {
                return node;
            }
        }
        return null;
    }

    @Nullable
    private static String determinePlatinumVersion(@NotNull OrangeConfig config) throws IOException {
        String pipelineVersionFile = config.pipelineVersionFile();
        if (pipelineVersionFile == null) {
            LOGGER.warn("No platinum version could be determined as pipeline version file was not passed");
            return null;
        }

        String platinumVersion = PipelineVersionFile.majorDotMinorVersion(pipelineVersionFile);
        if (platinumVersion != null) {
            LOGGER.info("Determined platinum version to be 'v{}'", platinumVersion);
        } else {
            LOGGER.warn("No platinum version could be determined as version could not be resolved from {}", pipelineVersionFile);
        }
        return platinumVersion;
    }

    @Nullable
    private static OrangeSample loadSampleData(@NotNull OrangeConfig config, boolean loadTumorSample) throws IOException {
        if (loadTumorSample) {
            LOGGER.info("Loading tumor sample data");
        } else {
            if (config.refSampleWGSMetricsFile() != null && config.refSampleFlagstatFile() != null) {
                LOGGER.info("Loading reference sample data");
            } else {
                LOGGER.info("Skipping loading of reference sample data as no flagstat or WGS metrics has been provided");
                return null;
            }
        }

        String metricsFile = loadTumorSample ? config.tumorSampleWGSMetricsFile() : config.refSampleWGSMetricsFile();
        WGSMetrics metrics = WGSMetricsFile.read(metricsFile);
        LOGGER.info(" Loaded WGS metrics from {}", metricsFile);

        String flagstatFile = loadTumorSample ? config.tumorSampleFlagstatFile() : config.refSampleFlagstatFile();
        Flagstat flagstat = FlagstatFile.read(flagstatFile);
        LOGGER.info(" Loaded flagstat from {}", flagstatFile);

        return ImmutableOrangeSample.builder().metrics(metrics).flagstat(flagstat).build();
    }

    @Nullable
    private static Map<String, Double> loadGermlineMVLHPerGene(@NotNull OrangeConfig config) throws IOException {
        String sageGermlineGeneCoverageTsv = config.sageGermlineGeneCoverageTsv();
        if (sageGermlineGeneCoverageTsv == null) {
            LOGGER.info("Skipping loading of germline MVLH as no germline gene coverage has been provided");
            return null;
        }

        Map<String, Double> mvlhPerGene = Maps.newTreeMap();

        List<String> lines = Files.readAllLines(new File(sageGermlineGeneCoverageTsv).toPath());
        String header = lines.get(0);

        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, GeneDepthFile.DELIM);
        int geneIndex = fieldsIndexMap.get(GeneDepthFile.COL_GENE);
        int mvlhIndex = fieldsIndexMap.get(GeneDepthFile.COL_MV_LIKELIHOOD);

        for (String line : lines.subList(1, lines.size())) {
            String[] values = line.split(GeneDepthFile.DELIM);
            String gene = values[geneIndex];
            String mvlhString = values[mvlhIndex].substring(0, values[mvlhIndex].length() - 1);
            double missedVariantLikelihood = Double.parseDouble(mvlhString) / 100D;
            mvlhPerGene.put(gene, missedVariantLikelihood);
        }
        return mvlhPerGene;
    }

    @NotNull
    private static PurpleData loadPurpleData(@NotNull OrangeConfig config) throws IOException {
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
        LOGGER.info(" Loaded {} gene copy numbers entries", purple.allSomaticGeneCopyNumbers().size());

        if (referenceSample != null) {
            LOGGER.info(" Loaded {} germline driver catalog entries", purple.germlineDrivers().size());
            LOGGER.info(" Loaded {} germline variants (of which {} are reportable)",
                    purple.allGermlineVariants().size(),
                    purple.reportableGermlineVariants().size());

            LOGGER.info(" Loaded {} germline deletions (of which {} are reportable)",
                    purple.allGermlineDeletions().size(),
                    purple.reportableGermlineDeletions().size());
        } else {
            LOGGER.debug(" Skipped loading germline variants since no reference sample configured");
        }

        return purple;
    }

    @NotNull
    private static LinxData loadLinxData(@NotNull OrangeConfig config) throws IOException {
        String linxGermlineDataDirectory = config.linxGermlineDataDirectory();
        LOGGER.info("Loading LINX somatic data from {}", config.linxSomaticDataDirectory());

        LinxData linx = LinxDataLoader.load(config.tumorSampleId(), config.linxSomaticDataDirectory(), linxGermlineDataDirectory);

        LOGGER.info(" Loaded {} structural variants", linx.allStructuralVariants().size());
        LOGGER.info(" Loaded {} structural drivers", linx.drivers().size());
        LOGGER.info(" Loaded {} fusions (of which {} are reportable)", linx.allFusions().size(), linx.reportableFusions().size());
        LOGGER.info(" Loaded {} breakends (of which {} are reportable)", linx.allBreakends().size(), linx.reportableBreakends().size());
        LOGGER.info(" Loaded {} reportable homozygous disruptions", linx.homozygousDisruptions().size());

        if (linxGermlineDataDirectory != null) {
            LOGGER.info("Loading LINX germline data from {}", linxGermlineDataDirectory);
            LOGGER.info(" Loaded {} germline disruptions (of which {} are reportable)",
                    linx.allGermlineDisruptions().size(),
                    linx.reportableGermlineDisruptions().size());
        } else {
            LOGGER.info(" Skipped loading LINX germline data as no linx germline data directory has been provided");
        }

        return linx;
    }

    @Nullable
    private IsofoxInterpretedData loadIsofoxData(@NotNull OrangeConfig config, @NotNull LinxInterpretedData linx) throws IOException {
        OrangeRNAConfig rna = config.rnaConfig();
        if (rna == null) {
            LOGGER.info("Skipping ISOFOX data loading as RNA is not configured");
            return null;
        }

        String isofoxCancerType = cohortMapper.cancerTypeForSample(createSample(config));
        if (isofoxCancerType == null) {
            LOGGER.warn("Could not resolve isofox cancer type for {}" + config.tumorSampleId());
            return null;
        }

        return IsofoxInterpreter.interpret(IsofoxDataLoader.load(isofoxCancerType,
                rna.isofoxGeneDistributionCsv(),
                rna.isofoxAltSjCohortCsv(),
                rna.isofoxSummaryCsv(),
                rna.isofoxGeneDataCsv(),
                rna.isofoxFusionCsv(),
                rna.isofoxAltSpliceJunctionCsv()), linx, driverGenes, knownFusionCache);
    }

    @NotNull
    private static LilacSummaryData loadLilacData(@NotNull OrangeConfig config) throws IOException {
        return LilacSummaryData.load(config.lilacQcCsv(), config.lilacResultCsv());
    }

    @NotNull
    private static VirusInterpreterData loadVirusInterpreterData(@NotNull OrangeConfig config) throws IOException {
        return VirusInterpreterDataLoader.load(config.annotatedVirusTsv());
    }

    @NotNull
    private static ChordData loadChordAnalysis(@NotNull OrangeConfig config) throws IOException {
        LOGGER.info("Loading CHORD data from {}", new File(config.chordPredictionTxt()).getParent());
        ChordData chordData = ChordDataFile.read(config.chordPredictionTxt());
        LOGGER.info(" HR Status: {} with type '{}'", chordData.hrStatus().display(), chordData.hrdType());
        return chordData;
    }

    @NotNull
    private static CuppaData loadCuppaData(@NotNull OrangeConfig config) throws IOException {
        LOGGER.info("Loading CUPPA from {}", new File(config.cuppaResultCsv()).getParent());
        List<CuppaDataFile> cuppaEntries = CuppaDataFile.read(config.cuppaResultCsv());
        LOGGER.info(" Loaded {} entries from {}", cuppaEntries.size(), config.cuppaResultCsv());

        CuppaData cuppaData = CuppaDataFactory.create(cuppaEntries);
        CuppaPrediction best = cuppaData.predictions().get(0);
        LOGGER.info(" Predicted cancer type '{}' with likelihood {}", best.cancerType(), best.likelihood());

        return cuppaData;
    }

    @Nullable
    private static List<PeachGenotype> loadPeachData(@NotNull OrangeConfig config) throws IOException {
        String peachGenotypeTsv = config.peachGenotypeTsv();

        if (peachGenotypeTsv == null) {
            LOGGER.info("Skipping PEACH loading since no peach genotype tsv has been provided");
            return null;
        }

        LOGGER.info("Loading PEACH from {}", new File(peachGenotypeTsv).getParent());
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(config.peachGenotypeTsv());
        LOGGER.info(" Loaded {} PEACH genotypes from {}", peachGenotypes.size(), config.peachGenotypeTsv());

        return peachGenotypes;
    }

    @NotNull
    private Map<PercentileType, Evaluation> evaluateCohortPercentiles(@NotNull OrangeConfig config, @NotNull PurpleInterpretedData purple) {
        PercentileType type = PercentileType.SV_TMB;

        Observation svTmbObservation = ImmutableObservation.builder()
                .sample(createSample(config))
                .type(type)
                .value(purple.characteristics().svTumorMutationalBurden())
                .build();

        LOGGER.info("Determining SV TMB percentile for value {}", svTmbObservation.value());
        Map<PercentileType, Evaluation> evaluations = Maps.newHashMap();
        Evaluation evaluation = percentilesModel.percentile(svTmbObservation);
        if (evaluation != null) {
            String cancerType = evaluation.cancerType();
            Double cancerTypePercentile = evaluation.cancerTypePercentile();
            LOGGER.info(" Determined percentile '{}' for pan-cancer and '{}' for cancer type '{}'",
                    evaluation.panCancerPercentile(),
                    cancerTypePercentile != null ? cancerTypePercentile : "NA",
                    cancerType != null ? cancerType : "NA");
            evaluations.put(type, evaluation);
        } else {
            LOGGER.warn("Could not evaluate SV TMB percentile for {}!", config.tumorSampleId());
        }

        return evaluations;
    }

    @NotNull
    private static OrangePlots buildPlots(@NotNull OrangeConfig config) {
        LOGGER.info("Loading plots");
        String linxPlotDir = config.linxPlotDirectory();
        List<String> linxDriverPlots = Lists.newArrayList();
        if (new File(linxPlotDir).exists()) {
            for (String file : new File(linxPlotDir).list()) {
                linxDriverPlots.add(linxPlotDir + File.separator + file);
            }
            LOGGER.info(" Loaded {} linx plots from {}", linxDriverPlots.size(), linxPlotDir);
        } else {
            LOGGER.debug(" No linx plots have been loaded as plot directory {} does not exist", linxPlotDir);
        }

        String kataegisPlot = config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".somatic.rainfall.png";
        if (!new File(kataegisPlot).exists()) {
            LOGGER.debug(" Could not locate kataegis plot '{}'", kataegisPlot);
            kataegisPlot = null;
        }

        return ImmutableOrangePlots.builder()
                .sageReferenceBQRPlot(config.sageSomaticRefSampleBQRPlot())
                .sageTumorBQRPlot(config.sageSomaticTumorSampleBQRPlot())
                .purpleInputPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".input.png")
                .purpleFinalCircosPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".circos.png")
                .purpleClonalityPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".somatic.clonality.png")
                .purpleCopyNumberPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".copynumber.png")
                .purpleVariantCopyNumberPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".somatic.png")
                .purplePurityRangePlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".purity.range.png")
                .purpleKataegisPlot(kataegisPlot)
                .linxDriverPlots(linxDriverPlots)
                .cuppaSummaryPlot(config.cuppaSummaryPlot())
                .cuppaFeaturePlot(config.cuppaFeaturePlot())
                .build();
    }

    @NotNull
    private static Sample createSample(@NotNull OrangeConfig config) {
        return ImmutableSample.builder().sampleId(config.tumorSampleId()).doids(config.primaryTumorDoids()).build();
    }
}
