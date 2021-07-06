package com.hartwig.hmftools.orange.algo;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOriginFile;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.pipeline.PipelineVersionFile;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.protect.chord.ChordDataLoader;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.virusinterpreter.VirusInterpreterData;
import com.hartwig.hmftools.protect.virusinterpreter.VirusInterpreterDataLoader;

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
        OrangePlots plots = ImmutableOrangePlots.builder().purpleCircosPlot(config.purpleCircosPlot()).build();

        return ImmutableOrangeReport.builder()
                .sampleId(config.tumorSampleId())
                .pipelineVersion(loadPipelineVersion(config))
                .configuredPrimaryTumor(loadConfiguredPrimaryTumor(config))
                .cuppaPrimaryTumor(loadCuppaPrimaryTumor(config))
                .purple(loadPurpleData(config))
                .linx(loadLinxData(config))
                .virusInterpreter(loadVirusInterpreterData(config))
                .chord(loadChordAnalysis(config))
                .plots(plots)
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
}
