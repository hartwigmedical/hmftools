package com.hartwig.hmftools.protect.chord;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.protect.ProtectConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ChordDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(ChordDataLoader.class);

    private ChordDataLoader() {
    }

    @NotNull
    public static ChordAnalysis load(@NotNull ProtectConfig config) throws IOException {
        LOGGER.info("Loading CHORD data from {}", new File(config.chordPredictionTxt()).getParent());
        ChordAnalysis chordAnalysis = ChordFileReader.read(config.chordPredictionTxt());
        LOGGER.info(" HR Status: {}", chordAnalysis.hrStatus().display());
        return chordAnalysis;
    }
}
