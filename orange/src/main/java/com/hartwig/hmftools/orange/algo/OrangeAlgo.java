package com.hartwig.hmftools.orange.algo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordDataLoader;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOriginFile;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDataLoader;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleDataLoader;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.orange.OrangeConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class OrangeAlgo {

    private static final Logger LOGGER = LogManager.getLogger(OrangeAlgo.class);

    @NotNull
    private final DoidEntry doidEntry;

    @NotNull
    public static OrangeAlgo fromConfig(@NotNull OrangeConfig config) throws IOException {
        LOGGER.info("Loading DOID database from {}", config.doidJsonFile());
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJsonFile());
        return new OrangeAlgo(doidEntry);
    }

    private OrangeAlgo(@NotNull final DoidEntry doidEntry) {
        this.doidEntry = doidEntry;
    }

    @NotNull
    public OrangeReport run(@NotNull OrangeConfig config) throws IOException {
        return ImmutableOrangeReport.builder()
                .sampleId(config.tumorSampleId())
                .pipelineVersion(loadPipelineVersion(config))
                .configuredPrimaryTumor(loadConfiguredPrimaryTumor(config))
                .cuppaPrimaryTumor(loadCuppaPrimaryTumor(config))
                .purple(loadPurpleData(config))
                .linx(loadLinxData(config))
                .virusInterpreter(loadVirusInterpreterData(config))
                .chord(loadChordAnalysis(config))
                .protect(loadProtectData(config))
                .germlineMVLHPerGene(loadGermlineMVLHPerGene(config))
                .plots(buildPlots(config))
                .build();
    }

    @NotNull
    private static String loadPipelineVersion(@NotNull OrangeConfig config) throws IOException {
        String pipelineVersion = PipelineVersionFile.majorDotMinorVersion(config.pipelineVersionFile());
        LOGGER.info("Loaded pipeline version '{}'", pipelineVersion);
        return pipelineVersion;
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

    @NotNull
    private static String loadCuppaPrimaryTumor(@NotNull OrangeConfig config) throws IOException {
        LOGGER.info("Loading Cuppa result from {}", new File(config.cuppaConclusionTxt()).getParent());
        String cuppaTumorLocation = MolecularTissueOriginFile.read(config.cuppaConclusionTxt()).conclusion();
        LOGGER.info(" Cuppa predicted primary tumor: {}", cuppaTumorLocation);
        return cuppaTumorLocation;
    }

    @NotNull
    private static PurpleData loadPurpleData(@NotNull OrangeConfig config) throws IOException {
        return PurpleDataLoader.load(config.tumorSampleId(),
                config.referenceSampleId(),
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleSomaticDriverCatalogTsv(),
                config.purpleSomaticVariantVcf(),
                config.purpleGermlineDriverCatalogTsv(),
                config.purpleGermlineVariantVcf(),
                config.purpleGeneCopyNumberTsv(),
                null,
                RefGenomeVersion.V37);
    }

    @NotNull
    private static LinxData loadLinxData(@NotNull OrangeConfig config) throws IOException {
        return LinxDataLoader.load(config.linxFusionTsv(), config.linxBreakendTsv(), config.linxDriverCatalogTsv());
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
    private static List<ProtectEvidence> loadProtectData(@NotNull OrangeConfig config) throws IOException {
        return ProtectEvidenceFile.read(config.protectEvidenceTsv());
    }

    @NotNull
    private static Map<String, Double> loadGermlineMVLHPerGene(@NotNull OrangeConfig config) throws IOException {
        Map<String, Double> mvlhPerGene = Maps.newHashMap();
        List<String> lines = Files.readAllLines(new File(config.sageGermlineGeneCoverageTsv()).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] values = line.split("\t");
            double mvlh = Double.parseDouble(values[1].substring(0, values[1].length() -1)) / 100D;
            mvlhPerGene.put(values[0], mvlh);
        }
        return mvlhPerGene;
    }

    @NotNull
    private static OrangePlots buildPlots(@NotNull OrangeConfig config) {
        LOGGER.info("Loading plots");
        String linxPlotDir = config.linxPlotDirectory();
        List<String> linxDriverPlots = Lists.newArrayList();
        for (String file : new File(linxPlotDir).list()) {
            linxDriverPlots.add(linxPlotDir + File.separator + file);
        }
        LOGGER.info(" Loaded {} linx plots from {}", linxDriverPlots.size(), linxPlotDir);

        return ImmutableOrangePlots.builder()
                .sageReferenceBQRPlot(config.sageSomaticRefSampleBQRPlot())
                .sageTumorBQRPlot(config.sageSomaticTumorSampleBQRPlot())
                .purpleInputPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".input.png")
                .purpleFinalCircosPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".circos.png")
                .purpleClonalityPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".somatic.clonality.png")
                .purpleCopyNumberPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".copynumber.png")
                .purpleVariantCopyNumberPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".somatic.png")
                .purplePurityRangePlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".purity.range.png")
                .purpleKataegisPlot(config.purplePlotDirectory() + File.separator + config.tumorSampleId() + ".somatic.rainfall.png")
                .linxDriverPlots(linxDriverPlots)
                .build();
    }
}
