package com.hartwig.hmftools.common.chord;

import java.io.File;
import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ChordDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(ChordDataLoader.class);

    private ChordDataLoader() {
    }

    @NotNull
    public static ChordAnalysis load(@NotNull String chordPredictionTxt) throws IOException {
        LOGGER.info("Loading CHORD data from {}", new File(chordPredictionTxt).getParent());
        ChordAnalysis chordAnalysis = ChordFileReader.read(chordPredictionTxt);
        LOGGER.info(" HR Status: {} with type '{}'", chordAnalysis.hrStatus().display(), chordAnalysis.hrdType());
        return chordAnalysis;
    }
}
