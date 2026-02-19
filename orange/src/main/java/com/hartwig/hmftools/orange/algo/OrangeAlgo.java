package com.hartwig.hmftools.orange.algo;

import static com.hartwig.hmftools.common.metrics.GeneDepthFile.generateGeneCoverageFilename;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.redux.BqrFile;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeDoidNode;
import com.hartwig.hmftools.datamodel.orange.OrangeDoidNode;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxData;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxDataLoader;
import com.hartwig.hmftools.orange.algo.linx.LinxData;
import com.hartwig.hmftools.orange.algo.linx.LinxDataLoader;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
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
import com.hartwig.hmftools.orange.algo.cuppa.CuppaDataFactory;
import com.hartwig.hmftools.orange.algo.immuno.ImmuneEscapeInterpreter;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpreter;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpreter;
import com.hartwig.hmftools.orange.algo.linx.LinxReportableClusters;
import com.hartwig.hmftools.orange.algo.plot.DummyPlotManager;
import com.hartwig.hmftools.orange.algo.plot.FileBasedPlotManager;
import com.hartwig.hmftools.orange.algo.plot.PlotManager;
import com.hartwig.hmftools.orange.algo.purple.PurpleData;
import com.hartwig.hmftools.orange.algo.purple.PurpleDataLoader;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpreter;
import com.hartwig.hmftools.orange.algo.sage.GermlineMVLHFactory;
import com.hartwig.hmftools.orange.algo.sigs.SigsEtiologiesLoader;
import com.hartwig.hmftools.orange.algo.sigs.SigsInterpreter;
import com.hartwig.hmftools.orange.algo.util.ReportLimiter;
import com.hartwig.hmftools.orange.algo.virus.VirusInterpreter;
import com.hartwig.hmftools.orange.algo.wildtype.WildTypeAlgo;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableSample;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.OrangeConversion;

import org.jetbrains.annotations.Nullable;

public class OrangeAlgo
{
    private final Map<String,DriverGene> mDriverGenes;
    private final Map<String, String> mEtiologyPerSignature;
    private final PlotManager mPlotManager;

    public static OrangeAlgo fromConfig(final OrangeConfig config) throws IOException
    {
        LOGGER.info("Reading driver genes from {}", config.DriverGenePanelTsv);
        List<DriverGene> driverGenes = DriverGeneFile.read(config.DriverGenePanelTsv);
        LOGGER.info(" Read {} driver genes", driverGenes.size());

        LOGGER.info("Reading signatures etiology from {}", config.SignaturesEtiologyTsv);
        Map<String, String> etiologyPerSignature = SigsEtiologiesLoader.read(config.SignaturesEtiologyTsv);
        LOGGER.info(" Read {} signatures etiology", etiologyPerSignature.size());

        String outputDir = config.OutputDir;
        PlotManager plotManager = !outputDir.isEmpty() ? new FileBasedPlotManager(outputDir) : new DummyPlotManager();

        return new OrangeAlgo(driverGenes, etiologyPerSignature, plotManager);
    }

    private OrangeAlgo(
            final List<DriverGene> driverGenes, final Map<String, String> etiologyPerSignature, final PlotManager plotManager)
    {
        mDriverGenes = Maps.newHashMap();
        driverGenes.forEach(x -> mDriverGenes.put(x.gene(), x));

        mEtiologyPerSignature = etiologyPerSignature;
        mPlotManager = plotManager;
    }

    public OrangeRecord run(final OrangeConfig config) throws Exception
    {
        Set<OrangeDoidNode> primaryTumorDoids = formConfiguredPrimaryTumorDoid(config);

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

        CytoBands cytoBands = new CytoBands(config.RefGenVersion);

        LinxInterpreter linxInterpreter = new LinxInterpreter(cytoBands);

        LinxRecord linx = linxInterpreter.interpret(linxData);

        List<DriverGene> driverGenes = new ArrayList<>(mDriverGenes.values());

        PurpleInterpreter purpleInterpreter = new PurpleInterpreter();

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

        boolean hasRefSample = config.ReferenceId != null;

        OrangeRecord orangeRecord = ImmutableOrangeRecord.builder()
                .sampleId(config.TumorId)
                .samplingDate(config.SamplingDate)
                .experimentType(config.RunType)
                .configuredPrimaryTumor(primaryTumorDoids)
                .refGenomeVersion(config.orangeRefGenomeVersion())
                .pipelineVersion(pipelineVersion)
                .refSample(refSample)
                .tumorSample(tumorSample)
                .germlineMVLHPerGene(mvlhPerGene)
                .purple(purple)
                .linx(linx)
                .wildTypeGenes(wildTypeGenes)
                .isofox(isofox)
                .lilac(lilac != null ? OrangeConversion.convert(lilac, hasRefSample, config.RnaSampleId != null) : null)
                .immuneEscape(immuneEscape)
                .virusInterpreter(virusInterpreter != null ? VirusInterpreter.interpret(virusInterpreter) : null)
                .chord(chord != null ? OrangeConversion.convert(chord) : null)
                .cuppa(cuppa)
                .peach(ConversionUtil.mapToIterable(peach, OrangeConversion::convert))
                .sigAllocations(SigsInterpreter.interpret(sigAllocations, mEtiologyPerSignature))
                .plots(buildPlots(config))
                .build();

        verifyPlots(orangeRecord.plots(), linxData);

        if(config.LimitJsonOutput)
        {
            orangeRecord = ReportLimiter.limitAllListsToMaxOne(orangeRecord);
        }

        return orangeRecord;
    }

    private Set<OrangeDoidNode> formConfiguredPrimaryTumorDoid(final OrangeConfig config)
    {
        Set<OrangeDoidNode> orangeNodes = Sets.newHashSet();

        if(!config.PrimaryTumorDoids.isEmpty())
        {
            /* currently unused
            DoidEntry doidEntry = null;

            if(config.DoidJsonFile != null)
            {
                LOGGER.info("Loading DOID database from {}", config.DoidJsonFile);
                doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(config.DoidJsonFile);
            }
            */

            LOGGER.debug("Determining configured primary tumor");

            Set<DoidNode> nodes = Sets.newHashSet();
            for(String doid : config.PrimaryTumorDoids)
            {
                DoidNode matchedNode = null;

                for(DoidNode node : nodes)
                {
                    if(node.doid().equals(doid))
                    {
                        matchedNode = node;
                        break;
                    }
                }

                if(matchedNode != null)
                {
                    LOGGER.debug(" Adding DOID {} ({}) as configured primary tumor", doid, matchedNode.doidTerm());
                    orangeNodes.add(OrangeConversion.convert(matchedNode));
                }
                else
                {
                    LOGGER.warn("Could not resolve doid '{}'", doid);
                }
            }
        }
        else if(config.PrimaryTumorLocation != null)
        {
            OrangeDoidNode doidNode = ImmutableOrangeDoidNode.builder().doid("").doidTerm(config.PrimaryTumorLocation).build();
            orangeNodes.add(doidNode);
        }

        return orangeNodes;
    }

    @Nullable
    private static String determinePipelineVersion(final OrangeConfig config) throws IOException
    {
        String pipelineVersionFile = config.PipelineVersionFile;
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
            LOGGER.info("Loading tumor metrics data");
        }
        else
        {
            if(config.ReferenceId != null)
            {
                LOGGER.info("Loading reference metrics data");
            }
            else
            {
                return null;
            }
        }

        String metricsFile = loadTumorSample ?
                BamMetricSummary.generateFilename(config.TumorBamMetricsDir, config.TumorId) :
                BamMetricSummary.generateFilename(config.ReferenceBamMetricsDir, config.ReferenceId);

        String flagstatFile = loadTumorSample ?
                BamFlagStats.generateFilename(config.TumorBamMetricsDir, config.TumorId) :
                BamFlagStats.generateFilename(config.ReferenceBamMetricsDir, config.ReferenceId);

        WGSMetrics metrics = OrangeConversion.convert(BamMetricSummary.read(metricsFile));
        LOGGER.info(" Loaded bam metrics from {}", metricsFile);

        Flagstat flagstat = OrangeConversion.convert(BamFlagStats.read(flagstatFile));
        LOGGER.info(" Loaded bam flagstats from {}", flagstatFile);

        return ImmutableOrangeSample.builder().metrics(metrics).flagstat(flagstat).build();
    }

    @Nullable
    private static Map<String, Double> loadGermlineMVLHPerGene(final OrangeConfig config, final Map<String,DriverGene> driverGenes)
            throws IOException
    {
        if(config.ReferenceId == null)
            return null;

        String geneCoverageFile = generateGeneCoverageFilename(config.ReferenceBamMetricsDir, config.ReferenceId);
        if(!Files.exists(Paths.get(geneCoverageFile)))
        {
            LOGGER.info("Skipping loading of germline MVLH as no germline gene coverage has been provided");
            return null;
        }

        Map<String, Double> mvlhPerGene = GermlineMVLHFactory.loadGermlineMVLHPerGene(geneCoverageFile, driverGenes);
        LOGGER.info("Loaded MVLH data for {} genes", mvlhPerGene.keySet().size());

        return mvlhPerGene;
    }

    private PurpleData loadPurpleData(final OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading PURPLE data from {}", config.PurpleDataDirectory);

        PurpleData purple = PurpleDataLoader.load(config, mDriverGenes);
        LOGGER.info(" Loaded {} somatic driver catalog entries", purple.somaticDrivers().size());
        LOGGER.info(" Loaded {} somatic variants", purple.somaticVariants().size());
        LOGGER.info(" Loaded {} somatic copy numbers entries", purple.somaticCopyNumbers().size());
        LOGGER.info(" Loaded {} somatic gene copy numbers entries", purple.somaticGeneCopyNumbers().size());

        if(config.ReferenceId != null)
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
        LOGGER.info("Loading LINX somatic data from {}", config.LinxSomaticDataDirectory);

        String linxGermlineDataDirectory = config.ReferenceId != null ? config.LinxGermlineDataDirectory : null;

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
        if(config.RnaSampleId == null || config.IsofoxDir == null)
        {
            LOGGER.info("Skipping Isofox data loading as RNA is not configured");
            return null;
        }

        return IsofoxDataLoader.load(config.TumorId, config.IsofoxDir);
    }

    @Nullable
    private static LilacSummaryData loadLilacData(final OrangeConfig config) throws IOException
    {
        if(config.LilacDir == null || !Files.exists(Paths.get(config.LilacDir)))
        {
            LOGGER.info("Skipping loading Lilac results since input dir not provided");
            return null;
        }

        LOGGER.info("Loading Lilac data from {}", config.LilacDir);

        return LilacSummaryData.read(config.LilacDir, config.TumorId);
    }

    @Nullable
    private static VirusInterpreterData loadVirusInterpreterData(final OrangeConfig config) throws IOException
    {
        if(config.VirusDir == null || !Files.exists(Paths.get(config.VirusDir)))
        {
            LOGGER.debug("Skipping Virus annotations as no input has been provided");
            return null;
        }

        return VirusInterpreterDataLoader.load(config.VirusDir, config.TumorId);
    }

    @Nullable
    private static ChordData loadChordAnalysis(final OrangeConfig config) throws IOException
    {
        if(config.ChordDir == null || !Files.exists(Paths.get(config.ChordDir)))
        {
            LOGGER.debug("Skipping Chord loading as no input has been provided");
            return null;
        }

        String chordFile = ChordDataFile.generateFilename(config.ChordDir, config.TumorId);
        LOGGER.info("Loading CHORD data from {}", chordFile);
        ChordData chordData = ChordDataFile.read(chordFile);
        LOGGER.info(" HR Status: {} with type '{}'", chordData.hrStatus().display(), chordData.hrdType());
        return chordData;
    }

    @Nullable
    private static CuppaData loadCuppaData(final OrangeConfig config) throws Exception
    {
        if(config.ReferenceId == null) // may change, esp if vCuppa is used
            return null;

        if(config.CuppaDir == null || !Files.exists(Paths.get(config.CuppaDir)))
        {
            LOGGER.debug("Skipping Cuppa loading as no input has been provided");
            return null;
        }

        String cuppaVisDataTsv = CuppaPredictions.generateVisDataTsvFilename(config.CuppaDir, config.TumorId);
        LOGGER.info("Loading CUPPA predictions from {}", cuppaVisDataTsv);
        CuppaData cuppaData = CuppaDataFactory.create(cuppaVisDataTsv);
        LOGGER.info(" Loaded {} CUPPA predictions from {}", cuppaData.predictions().size(), cuppaVisDataTsv);

        return cuppaData;
    }

    @Nullable
    private static List<PeachGenotype> loadPeachData(final OrangeConfig config) throws IOException
    {
        if(config.ReferenceId == null)
            return null;

        if(config.PeachDir == null || !Files.exists(Paths.get(config.PeachDir)))
        {
            LOGGER.debug("Skipping Peach loading as no input has been provided");
            return null;
        }

        String peachGenotypeFile = PeachGenotypeFile.generateFileName(config.PeachDir, config.ReferenceId);

        LOGGER.info("Loading Peach from {}", peachGenotypeFile);
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(peachGenotypeFile);
        LOGGER.info(" Loaded {} Peach genotypes from {}", peachGenotypes.size(), peachGenotypeFile);
        List<PeachGenotype> filterUGT1A1FromPeachGenotypes = peachGenotypes.stream().filter(genotype -> !genotype.gene().equals("UGT1A1")).toList();
        return filterUGT1A1FromPeachGenotypes;
    }

    @Nullable
    private static List<SignatureAllocation> loadSigAllocations(final OrangeConfig config) throws IOException
    {
        if(config.SigsDir == null || !Files.exists(Paths.get(config.SigsDir)))
        {
            LOGGER.debug("Skipping Signature loading as no input has been provided");
            return null;
        }

        String sigsAllocationFile = SignatureAllocationFile.generateFilename(config.SigsDir, config.TumorId);

        LOGGER.info("Loading Sigs from {}", sigsAllocationFile);
        List<SignatureAllocation> sigsAllocations = SignatureAllocationFile.read(sigsAllocationFile);
        LOGGER.info(" Loaded {} signature allocations from {}", sigsAllocations.size(), sigsAllocationFile);

        return sigsAllocations;
    }

    private OrangePlots buildPlots(final OrangeConfig config) throws IOException
    {
        LOGGER.info("Loading plots");

        mPlotManager.createPlotDirectory();

        String linxPlotDir = config.LinxPlotDirectory;
        List<String> linxDriverPlots = Lists.newArrayList();

        if(linxPlotDir != null)
        {
            for(String file : new File(linxPlotDir).list())
            {
                linxDriverPlots.add(mPlotManager.processPlotFile(linxPlotDir + File.separator + file));
            }

            LOGGER.info(" Loaded {} linx plots from {}", linxDriverPlots.size(), linxPlotDir);
        }

        String tumorBqrPlot = mPlotManager.processPlotFile(BqrFile.generatePlotFilename(config.TumorReduxDir, config.TumorId));

        if(tumorBqrPlot == null)
            tumorBqrPlot = "";

        String refBqrPlot = config.ReferenceId != null ?
                mPlotManager.processPlotFile(BqrFile.generatePlotFilename(config.ReferenceReduxDir, config.ReferenceId)) : "";

        String purplePlotBasePath = config.PurplePlotDirectory + File.separator + config.TumorId;
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
                config.ReferenceId != null ? CuppaPredictions.generateVisPlotFilename(config.CuppaDir, config.TumorId) : null);

        return ImmutableOrangePlots.builder()
                .sageReferenceBQRPlot(refBqrPlot)
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
                .sampleId(config.TumorId)
                .doids(config.PrimaryTumorDoids)
                .build();
    }
}
