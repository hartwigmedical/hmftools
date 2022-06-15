package com.hartwig.hmftools.orange.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.time.LocalDate;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordDataLoader;
import com.hartwig.hmftools.common.cuppa.CuppaData;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.CuppaEntry;
import com.hartwig.hmftools.common.cuppa.CuppaFactory;
import com.hartwig.hmftools.common.cuppa.CuppaPrediction;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.isofox.IsofoxDataLoader;
import com.hartwig.hmftools.common.lilac.LilacData;
import com.hartwig.hmftools.common.lilac.LilacDataLoader;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDataLoader;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleDataLoader;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeRNAConfig;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpretedData;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpreter;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpretedData;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpreter;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpreter;
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

    private static final DecimalFormat PERCENTAGE = new DecimalFormat("#'%'", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

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

        List<ProtectEvidence> protect = loadProtectData(config);
        LinxInterpretedData linx = LinxInterpreter.interpret(loadLinxData(config), protect, driverGenes);

        ChordAnalysis chord = loadChordAnalysis(config);
        PurpleInterpretedData purple = PurpleInterpreter.interpret(loadPurpleData(config), protect, driverGenes, chord);

        return ImmutableOrangeReport.builder()
                .sampleId(config.tumorSampleId())
                .reportDate(LocalDate.now())
                .configuredPrimaryTumor(configuredPrimaryTumor)
                .refGenomeVersion(config.refGenomeVersion())
                .platinumVersion(platinumVersion)
                .refSample(refSample)
                .tumorSample(tumorSample)
                .germlineMVLHPerGene(loadGermlineMVLHPerGene(config))
                .purple(purple)
                .linx(linx)
                .isofox(loadIsofoxData(config, linx))
                .lilac(loadLilacData(config))
                .virusInterpreter(loadVirusInterpreterData(config))
                .chord(chord)
                .cuppa(loadCuppaData(config))
                .peach(loadPeachData(config))
                .protect(protect)
                .cohortEvaluations(evaluateCohortPercentiles(config, purple))
                .plots(buildPlots(config))
                .build();
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
        if (pipelineVersionFile != null) {
            String platinumVersion = PipelineVersionFile.majorDotMinorVersion(pipelineVersionFile);
            if (platinumVersion != null) {
                LOGGER.info("Determined platinum version to be 'v{}'", platinumVersion);
            } else {
                LOGGER.warn("No platinum version could be determined as version could not be resolved from {}", pipelineVersionFile);
            }
            return platinumVersion;
        } else {
            LOGGER.warn("No platinum version could be determined as pipeline version file was not passed");
            return null;
        }
    }

    @NotNull
    private static OrangeSample loadSampleData(@NotNull OrangeConfig config, boolean loadTumorSample) throws IOException {
        if (loadTumorSample) {
            LOGGER.info("Loading tumor sample data");
        } else {
            LOGGER.info("Loading reference sample data");
        }

        String metricsFile = loadTumorSample ? config.tumorSampleWGSMetricsFile() : config.refSampleWGSMetricsFile();
        WGSMetrics metrics = WGSMetricsFile.read(metricsFile);
        LOGGER.info(" Loaded WGS metrics from {}", metricsFile);

        String flagstatFile = loadTumorSample ? config.tumorSampleFlagstatFile() : config.refSampleFlagstatFile();
        Flagstat flagstat = FlagstatFile.read(flagstatFile);
        LOGGER.info(" Loaded flagstat from {}", flagstatFile);

        return ImmutableOrangeSample.builder().metrics(metrics).flagstat(flagstat).build();
    }

    @NotNull
    private static Map<String, Double> loadGermlineMVLHPerGene(@NotNull OrangeConfig config) throws IOException {
        Map<String, Double> mvlhPerGene = Maps.newTreeMap();
        List<String> lines = Files.readAllLines(new File(config.sageGermlineGeneCoverageTsv()).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] values = line.split("\t");
            double mvlh = Double.parseDouble(values[1].substring(0, values[1].length() - 1)) / 100D;
            mvlhPerGene.put(values[0], mvlh);
        }
        return mvlhPerGene;
    }

    @NotNull
    private static PurpleData loadPurpleData(@NotNull OrangeConfig config) throws IOException {
        return PurpleDataLoader.load(config.tumorSampleId(),
                config.referenceSampleId(),
                config.rnaConfig() != null ? config.rnaConfig().rnaSampleId() : null,
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleSomaticDriverCatalogTsv(),
                config.purpleSomaticVariantVcf(),
                config.purpleGermlineDriverCatalogTsv(),
                config.purpleGermlineVariantVcf(),
                config.purpleGeneCopyNumberTsv(),
                config.purpleSomaticCopyNumberTsv(),
                config.purpleGermlineDeletionTsv(),
                config.refGenomeVersion());
    }

    @NotNull
    private static LinxData loadLinxData(@NotNull OrangeConfig config) throws IOException {
        return LinxDataLoader.load(config.linxStructuralVariantTsv(),
                config.linxFusionTsv(),
                config.linxBreakendTsv(),
                config.linxDriverCatalogTsv(),
                config.linxDriverTsv(),
                config.linxGermlineDisruptionTsv());
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
    private static LilacData loadLilacData(@NotNull OrangeConfig config) throws IOException {
        return LilacDataLoader.load(config.lilacQcCsv(), config.lilacResultCsv());
    }

    @NotNull
    private static VirusInterpreterData loadVirusInterpreterData(@NotNull OrangeConfig config) throws IOException {
        return VirusInterpreterDataLoader.load(config.annotatedVirusTsv());
    }

    @NotNull
    private static ChordAnalysis loadChordAnalysis(@NotNull OrangeConfig config) throws IOException {
        return ChordDataLoader.load(config.chordPredictionTxt());
    }

    @NotNull
    private static CuppaData loadCuppaData(@NotNull OrangeConfig config) throws IOException {
        LOGGER.info("Loading CUPPA from {}", new File(config.cuppaResultCsv()).getParent());
        List<CuppaEntry> cuppaEntries = CuppaDataFile.read(config.cuppaResultCsv());
        LOGGER.info(" Loaded {} entries from {}", cuppaEntries.size(), config.cuppaResultCsv());

        CuppaData cuppaData = CuppaFactory.create(cuppaEntries);
        CuppaPrediction best = cuppaData.predictions().get(0);
        LOGGER.info(" Predicted cancer type '{}' with likelihood {}", best.cancerType(), best.likelihood());

        return cuppaData;
    }

    @NotNull
    private static List<PeachGenotype> loadPeachData(@NotNull OrangeConfig config) throws IOException {
        LOGGER.info("Loading PEACH from {}", new File(config.peachGenotypeTsv()).getParent());
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(config.peachGenotypeTsv());
        LOGGER.info(" Loaded {} PEACH genotypes from {}", peachGenotypes.size(), config.peachGenotypeTsv());

        return peachGenotypes;
    }

    @NotNull
    private static List<ProtectEvidence> loadProtectData(@NotNull OrangeConfig config) throws IOException {
        LOGGER.info("Loading PROTECT data from {}", new File(config.protectEvidenceTsv()).getParent());
        List<ProtectEvidence> evidences = ProtectEvidenceFile.read(config.protectEvidenceTsv());
        LOGGER.info(" Loaded {} PROTECT evidences from {}", evidences.size(), config.protectEvidenceTsv());

        return evidences;
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
                    PERCENTAGE.format(evaluation.panCancerPercentile() * 100),
                    cancerTypePercentile != null ? PERCENTAGE.format(cancerTypePercentile * 100) : "NA",
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
