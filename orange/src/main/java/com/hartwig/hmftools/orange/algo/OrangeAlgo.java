package com.hartwig.hmftools.orange.algo;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.DEFAULT_PIPELINE_OUTPUT;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.algo.purple.PurpleDataLoader.addPurplePlots;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.common.sigs.SnvSigUtils;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeDoidNode;
import com.hartwig.hmftools.datamodel.orange.OrangeDoidNode;
import com.hartwig.hmftools.orange.algo.immuno.LilacInterpreter;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxData;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxDataLoader;
import com.hartwig.hmftools.orange.algo.linx.LinxData;
import com.hartwig.hmftools.orange.algo.linx.LinxDataLoader;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangePlots;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangePlots;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaDataFactory;
import com.hartwig.hmftools.orange.algo.immuno.ImmuneEscapeInterpreter;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpreter;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpreter;
import com.hartwig.hmftools.orange.algo.plot.DummyPlotManager;
import com.hartwig.hmftools.orange.algo.plot.FileBasedPlotManager;
import com.hartwig.hmftools.orange.algo.plot.PlotManager;
import com.hartwig.hmftools.orange.algo.purple.PurpleData;
import com.hartwig.hmftools.orange.algo.purple.PurpleDataLoader;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpreter;
import com.hartwig.hmftools.orange.algo.sigs.SigsInterpreter;
import com.hartwig.hmftools.orange.algo.virus.VirusInterpreter;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.OrangeConversion;

import org.jetbrains.annotations.Nullable;

public class OrangeAlgo
{
    private final Map<String,String> mEtiologyPerSignature;
    private final PlotManager mPlotManager;

    public static OrangeAlgo fromConfig(final OrangeConfig config) throws IOException
    {
        Map<String,String> etiologyPerSignature = SnvSigUtils.loadSnvSignatureEtiologies();

        String outputDir = config.OutputDir;
        PlotManager plotManager = !outputDir.isEmpty() ? new FileBasedPlotManager(outputDir) : new DummyPlotManager();

        return new OrangeAlgo(etiologyPerSignature, plotManager);
    }

    private OrangeAlgo(final Map<String,String> etiologyPerSignature, final PlotManager plotManager)
    {
        mEtiologyPerSignature = etiologyPerSignature;
        mPlotManager = plotManager;
    }

    public OrangeRecord run(final OrangeConfig config) throws Exception
    {
        LOGGER.info("building report data");

        Set<OrangeDoidNode> primaryTumorDoids = formConfiguredPrimaryTumorDoid(config);

        String pipelineVersion = determinePipelineVersion(config);

        PurpleData purpleData = loadPurpleData(config);
        LinxData linxData = loadLinxData(config);
        ChordData chord = loadChordAnalysis(config);
        LilacRecord lilac = loadLilacData(config);
        VirusInterpreterData virusInterpreter = loadVirusInterpreterData(config);
        CuppaData cuppa = loadCuppaData(config);
        List<PeachGenotype> peach = loadPeachData(config);
        List<SignatureAllocation> sigAllocations = loadSigAllocations(config);
        IsofoxData isofoxData = loadIsofoxData(config);

        CytoBands cytoBands = new CytoBands(config.RefGenVersion);

        LinxInterpreter linxInterpreter = new LinxInterpreter(cytoBands);

        LinxRecord linx = linxInterpreter.interpret(linxData, isofoxData);

        PurpleInterpreter purpleInterpreter = new PurpleInterpreter();

        PurpleRecord purple = purpleInterpreter.interpret(purpleData, linxData, isofoxData);

        ImmuneEscapeRecord immuneEscape = ImmuneEscapeInterpreter.interpret(purple, linx);

        IsofoxRecord isofox = null;
        if(isofoxData != null)
        {
            IsofoxInterpreter isofoxInterpreter = new IsofoxInterpreter(linx);
            isofox = isofoxInterpreter.interpret(isofoxData);
        }

        OrangePlots plots = buildPlots(config, linxData);

        OrangeRecord orangeRecord = ImmutableOrangeRecord.builder()
                .sampleId(config.TumorId)
                .referenceId(config.ReferenceId)
                .samplingDate(config.SamplingDate)
                .experimentType(config.RunType)
                .configuredPrimaryTumor(primaryTumorDoids)
                .refGenomeVersion(config.orangeRefGenomeVersion())
                .pipelineVersion(pipelineVersion)
                .purple(purple)
                .linx(linx)
                .isofox(isofox)
                .lilac(lilac)
                .immuneEscape(immuneEscape)
                .virusInterpreter(virusInterpreter != null ? VirusInterpreter.interpret(virusInterpreter) : null)
                .chord(chord != null ? OrangeConversion.convert(chord) : null)
                .cuppa(cuppa)
                .peach(ConversionUtil.mapToIterable(peach, OrangeConversion::convert))
                .sigAllocations(SigsInterpreter.interpret(sigAllocations, mEtiologyPerSignature))
                .plots(plots)
                .build();

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

            LOGGER.debug("determining configured primary tumor");

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

    private static String determinePipelineVersion(final OrangeConfig config) throws IOException
    {
        String defaultPipelineVersion = DEFAULT_PIPELINE_OUTPUT.toString();

        String pipelineVersionFile = config.PipelineVersionFile;
        if(pipelineVersionFile == null)
        {
            LOGGER.info("no pipeline version file, defaulting to {}", defaultPipelineVersion);
            return defaultPipelineVersion;
        }

        String pipelineVersion = PipelineVersionFile.majorDotMinorVersion(pipelineVersionFile);
        if(pipelineVersion != null)
        {
            LOGGER.info("determined pipeline version to be 'v{}'", pipelineVersion);
        }
        else
        {
            LOGGER.warn("no pipeline version determined from {}", pipelineVersionFile);
        }

        return defaultPipelineVersion;
    }

    private PurpleData loadPurpleData(final OrangeConfig config) throws IOException
    {
        PurpleData purple = PurpleDataLoader.load(config);
        LOGGER.debug(" loaded {} somatic driver catalog entries", purple.somaticDrivers().size());
        LOGGER.debug(" loaded {} somatic variants", purple.somaticVariants().size());

        if(config.ReferenceId != null)
        {
            LOGGER.debug(" loaded {} germline driver catalog entries", purple.germlineDrivers().size());
            LOGGER.debug(" loaded {} reportable germline variants", purple.germlineVariants().size());

            LOGGER.debug(" loaded {} reportable germline copy-number events", purple.germlineAmpDels().size());
        }
        else
        {
            LOGGER.debug(" Skipped loading germline variants and deletions since no reference sample configured");
        }

        return purple;
    }

    private static LinxData loadLinxData(final OrangeConfig config) throws IOException
    {
        String linxGermlineDataDirectory = config.ReferenceId != null ? config.LinxGermlineDataDirectory : null;

        LinxData linx = LinxDataLoader.load(config);

        LOGGER.debug(" loaded {} somatic structural variants", linx.somaticSvAnnotations().size());
        LOGGER.debug(" loaded {} somatic structural drivers", linx.somaticDriverData().size());
        LOGGER.debug(" loaded {} somatic fusions", linx.fusions().size());

        LOGGER.debug(" loaded {} somatic breakends", linx.somaticBreakends().size());

        if(linxGermlineDataDirectory != null)
        {
            LOGGER.debug(" loaded {} germline structural variants", linx.germlineSvAnnotations().size());
            LOGGER.debug(" loaded {} germline breakends (of which {} are reportable)",
                    linx.germlineBreakends().size(),
                    linx.germlineBreakends().size());
            LOGGER.debug(" loaded {} germline disruptions", linx.germlineDisruptions().size());
        }
        else if(config.hasReference())
        {
            LOGGER.warn(" skipped loading Linx germline data as no linx germline data directory has been provided");
        }

        return linx;
    }

    @Nullable
    private IsofoxData loadIsofoxData(final OrangeConfig config) throws IOException
    {
        if(config.RnaSampleId == null || config.IsofoxDir == null)
        {
            LOGGER.debug("skipping Isofox data loading as RNA is not configured");
            return null;
        }

        return IsofoxDataLoader.load(config.TumorId, config.IsofoxDir);
    }

    @Nullable
    private static LilacRecord loadLilacData(final OrangeConfig config) throws IOException
    {
        if(config.LilacDir == null || !Files.exists(Paths.get(config.LilacDir)))
        {
            LOGGER.info("skipping loading Lilac results since input dir not provided");
            return null;
        }

        LOGGER.debug("loading Lilac data from {}", config.LilacDir);

        return LilacInterpreter.build(config);
    }

    @Nullable
    private static VirusInterpreterData loadVirusInterpreterData(final OrangeConfig config) throws IOException
    {
        if(config.VirusDir == null || !Files.exists(Paths.get(config.VirusDir)))
        {
            LOGGER.debug("skipping Virus annotations as no input has been provided");
            return null;
        }

        return VirusInterpreterDataLoader.load(config.VirusDir, config.TumorId);
    }

    @Nullable
    private static ChordData loadChordAnalysis(final OrangeConfig config) throws IOException
    {
        if(!config.hasReference() || config.RunType == ExperimentType.TARGETED)
            return null;

        if(config.ChordDir == null || !Files.exists(Paths.get(config.ChordDir)))
        {
            LOGGER.debug("skipping Chord loading as no input has been provided");
            return null;
        }

        String chordFile = ChordDataFile.generateFilename(config.ChordDir, config.TumorId);
        ChordData chordData = ChordDataFile.read(chordFile);
        LOGGER.debug(" Chord HR Status: {} with type '{}'", chordData.hrStatus().display(), chordData.hrdType());
        return chordData;
    }

    @Nullable
    private static CuppaData loadCuppaData(final OrangeConfig config) throws Exception
    {
        if(config.ReferenceId == null || config.RunType == ExperimentType.TARGETED) // may change, esp if vCuppa is used
            return null;

        if(config.CuppaDir == null || !Files.exists(Paths.get(config.CuppaDir)))
        {
            LOGGER.debug("skipping Cuppa loading as no input has been provided");
            return null;
        }

        String cuppaVisDataTsv = CuppaPredictions.generateVisDataTsvFilename(config.CuppaDir, config.TumorId);
        CuppaData cuppaData = CuppaDataFactory.create(cuppaVisDataTsv);
        LOGGER.debug(" loaded {} Cuppa predictions from {}", cuppaData.predictions().size(), cuppaVisDataTsv);

        return cuppaData;
    }

    @Nullable
    private static List<PeachGenotype> loadPeachData(final OrangeConfig config) throws IOException
    {
        if(config.ReferenceId == null)
            return null;

        if(config.PeachDir == null || !Files.exists(Paths.get(config.PeachDir)))
        {
            LOGGER.debug("skipping Peach loading as no input has been provided");
            return null;
        }

        String peachGenotypeFile = PeachGenotypeFile.generateFileName(config.PeachDir, config.ReferenceId);

        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(peachGenotypeFile);
        LOGGER.debug(" loaded {} Peach genotypes from {}", peachGenotypes.size(), peachGenotypeFile);
        List<PeachGenotype> filterUGT1A1FromPeachGenotypes = peachGenotypes.stream().filter(genotype -> !genotype.gene().equals("UGT1A1")).toList();
        return filterUGT1A1FromPeachGenotypes;
    }

    @Nullable
    private static List<SignatureAllocation> loadSigAllocations(final OrangeConfig config) throws IOException
    {
        if(!config.hasReference() || config.RunType == ExperimentType.TARGETED)
            return null;

        if(config.SigsDir == null || !Files.exists(Paths.get(config.SigsDir)))
        {
            LOGGER.debug("skipping Signature loading as no input has been provided");
            return null;
        }

        String sigsAllocationFile = SignatureAllocationFile.generateFilename(config.SigsDir, config.TumorId);

        List<SignatureAllocation> sigsAllocations = SignatureAllocationFile.read(sigsAllocationFile);
        LOGGER.debug(" loaded {} signature allocations from {}", sigsAllocations.size(), sigsAllocationFile);

        return sigsAllocations;
    }

    private OrangePlots buildPlots(final OrangeConfig config, final LinxData linxData) throws IOException
    {
        LOGGER.debug("loading plots");

        mPlotManager.createPlotDirectory();

        ImmutableOrangePlots.Builder plotBuilder = ImmutableOrangePlots.builder();

        try
        {
            addPurplePlots(config, mPlotManager, plotBuilder);
        }
        catch(Exception e)
        {
            LOGGER.warn("failed to find required plots: {}", e.toString());
            System.exit(1);
        }

        List<String> linxDriverPlots = Lists.newArrayList();

        for(String linxPlot : linxData.reportableEventPlots())
        {
            linxDriverPlots.add(mPlotManager.processPlotFile(linxPlot));
        }

        plotBuilder.linxDriverPlots(linxDriverPlots);

        String qSeeSourcePlot = config.QSeeDirectory + File.separator + config.TumorId + ".qsee.vis.report.png";
        String qSeePlot = null;

        if(Files.exists(Paths.get(qSeeSourcePlot)))
        {
            qSeePlot = mPlotManager.processPlotFile(qSeeSourcePlot);
        }

        plotBuilder.qSeePlot(qSeePlot);

        String cuppaSummaryPlot = mPlotManager.processPlotFile(
                config.ReferenceId != null ? CuppaPredictions.generateVisPlotFilename(config.CuppaDir, config.TumorId) : null);

        plotBuilder.cuppaSummaryPlot(cuppaSummaryPlot);

        return plotBuilder.build();
    }
}
