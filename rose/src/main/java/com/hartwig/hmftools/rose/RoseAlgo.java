package com.hartwig.hmftools.rose;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordDataLoader;
import com.hartwig.hmftools.common.cuppa.ImmutableMolecularTissueOrigin;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOrigin;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOriginFile;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDataLoader;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleDataLoader;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.rose.actionability.ActionabilityEntry;
import com.hartwig.hmftools.rose.actionability.ActionabilityFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class RoseAlgo {
    private static final Logger LOGGER = LogManager.getLogger(RoseAlgo.class);

    @NotNull
    private final List<ActionabilityEntry> actionabilityEntry;
    @NotNull
    private final List<DriverGene> driverGenes;

    @NotNull
    public static RoseAlgo build(@NotNull String actionabilityDatabaseTsv, @NotNull String driverGene37Tsv,
            @NotNull String driverGene38Tsv, @NotNull RefGenomeVersion refGenomeVersion) throws IOException {
        List<ActionabilityEntry> actionabilityEntry = ActionabilityFileReader.read(actionabilityDatabaseTsv);
        List<DriverGene> driverGenes = refGenomeVersion == RefGenomeVersion.V37
                ? readDriverGenesFromFile(driverGene37Tsv)
                : readDriverGenesFromFile(driverGene38Tsv);
        return new RoseAlgo(actionabilityEntry, driverGenes);
    }

    private RoseAlgo(final @NotNull List<ActionabilityEntry> actionabilityEntry, final @NotNull List<DriverGene> driverGenes) {
        this.actionabilityEntry = actionabilityEntry;
        this.driverGenes = driverGenes;
    }

    @NotNull
    private static List<DriverGene> readDriverGenesFromFile(@NotNull String driverGeneTsv) throws IOException {
        LOGGER.info(" Reading driver genes from {}", driverGeneTsv);
        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsv);
        LOGGER.info("  Read {} driver gene entries", driverGenes.size());
        return driverGenes;
    }

    @NotNull
    public RoseData run(@NotNull RoseConfig config) throws IOException {
        PurpleData purple = loadPurpleData(config);

        return ImmutableRoseData.builder()
                .sampleId(config.tumorSampleId())
                .purple(purple)
                .linx(loadLinxData(config))
                .virusInterpreter(loadVirusInterpreterData(config))
                .chord(loadChordAnalysis(config))
                .molecularTissueOrigin(loadCuppaData(config))
                .actionabilityEntries(actionabilityEntry)
                .driverGenes(driverGenes)
                .build();
    }

    @NotNull
    private static PurpleData loadPurpleData(@NotNull RoseConfig config) throws IOException {
        return PurpleDataLoader.load(config.tumorSampleId(),
                config.refSampleId(),
                null,
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleSomaticDriverCatalogTsv(),
                config.purpleSomaticVariantVcf(),
                config.purpleGermlineDriverCatalogTsv(),
                config.purpleGermlineVariantVcf(),
                config.purpleSomaticCopyNumberTsv(),
                null,
                null,
                null);
    }

    @NotNull
    private static MolecularTissueOrigin loadCuppaData(@NotNull RoseConfig config) throws IOException {

        return ImmutableMolecularTissueOrigin.builder()
                .conclusion(MolecularTissueOriginFile.read(config.molecularTissueOriginTxt()).conclusion())
                .plotPath(Strings.EMPTY) //Plot isn't used in this tool
                .build();
    }

    @NotNull
    private static LinxData loadLinxData(@NotNull RoseConfig config) throws IOException {
        return LinxDataLoader.load(config.linxFusionTsv(), config.linxBreakendTsv(), null, config.linxDriverCatalogTsv(), null);
    }

    @NotNull
    private static VirusInterpreterData loadVirusInterpreterData(@NotNull RoseConfig config) throws IOException {
        return VirusInterpreterDataLoader.load(config.annotatedVirusTsv());
    }

    @NotNull
    private static ChordAnalysis loadChordAnalysis(@NotNull RoseConfig config) throws IOException {
        return ChordDataLoader.load(config.chordPredictionTxt());
    }
}