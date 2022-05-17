package com.hartwig.hmftools.summon;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordDataLoader;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDataLoader;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleDataLoader;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.summon.actionability.ActionabilityEntry;
import com.hartwig.hmftools.summon.actionability.ActionabilityFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SummonAlgo {
    private static final Logger LOGGER = LogManager.getLogger(SummonAlgo.class);

    @NotNull
    private final List<ActionabilityEntry> actionabilityEntry;

    @NotNull
    public static SummonAlgo build(@NotNull String actionabilityDatabaseTsv) throws IOException {
        List<ActionabilityEntry> actionabilityEntry = ActionabilityFileReader.read(actionabilityDatabaseTsv);

        return new SummonAlgo(actionabilityEntry);
    }

    private SummonAlgo(final @NotNull List<ActionabilityEntry> actionabilityEntry) {
        this.actionabilityEntry = actionabilityEntry;
    }

    @NotNull
    public SummonData run(@NotNull SummonConfig config) throws IOException {
        PurpleData purple = loadPurpleData(config);

        return ImmutableSummonData.builder()
                .sampleId(config.tumorSampleId())
                .purple(purple)
                .linx(loadLinxData(config))
                .virusInterpreter(loadVirusInterpreterData(config))
                .chord(loadChordAnalysis(config))
                .protect(loadProtectData(config))
                .actionabilityEntry(actionabilityEntry)
                .build();
    }

    @NotNull
    private static PurpleData loadPurpleData(@NotNull SummonConfig config) throws IOException {
        return PurpleDataLoader.load(config.tumorSampleId(),
                config.refSampleId(),
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleSomaticDriverCatalogTsv(),
                config.purpleSomaticVariantVcf(),
                config.purpleGermlineDriverCatalogTsv(),
                config.purpleGermlineVariantVcf(),
                config.purpleSomaticCopyNumberTsv());
    }

    @NotNull
    private static LinxData loadLinxData(@NotNull SummonConfig config) throws IOException {
        return LinxDataLoader.load(config.linxFusionTsv(),
                config.linxBreakendTsv(),
                null,
                config.linxDriverCatalogTsv(),
                null);
    }

    @NotNull
    private static VirusInterpreterData loadVirusInterpreterData(@NotNull SummonConfig config) throws IOException {
        return VirusInterpreterDataLoader.load(config.annotatedVirusTsv());
    }

    @NotNull
    private static ChordAnalysis loadChordAnalysis(@NotNull SummonConfig config) throws IOException {
        return ChordDataLoader.load(config.chordPredictionTxt());
    }

    @NotNull
    private static List<ProtectEvidence> loadProtectData(@NotNull SummonConfig config) throws IOException {
        LOGGER.info("Loading PROTECT data from {}", new File(config.protectEvidenceTsv()).getParent());
        List<ProtectEvidence> evidences = ProtectEvidenceFile.read(config.protectEvidenceTsv());
        LOGGER.info(" Loaded {} PROTECT evidences from {}", evidences.size(), config.protectEvidenceTsv());

        return evidences;
    }
}