package com.hartwig.hmftools.common.cuppa;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class MolecularTissueOriginFile {

    private static final Logger LOGGER = LogManager.getLogger(MolecularTissueOriginFile.class);

    public MolecularTissueOriginFile() {
    }

    @NotNull
    public static MolecularTissueOrginData read(@NotNull String molecularTissueOriginTxt) throws IOException {
        String origin = fromLines(Files.readAllLines(new File(molecularTissueOriginTxt).toPath()));
        if (origin == null) {
            LOGGER.warn("No molecular tissue origin could be read from {}!", molecularTissueOriginTxt);
            origin = Strings.EMPTY;
        }
        return extractPedictionDataOrigin(origin);
    }

    @Nullable
    private static String fromLines(@NotNull List<String> lines) {
        if (lines.size() == 1) {
            return lines.get(0).split(" - ")[1];
        } else {
            return null;
        }
    }

    @NotNull
    @VisibleForTesting
    public static MolecularTissueOrginData extractPedictionDataOrigin(@NotNull String cuppaResult) {
        String orgin = Strings.EMPTY;
        String prediction = null;
        String[] lengthCuppaResult = cuppaResult.split("\\(");

        if (lengthCuppaResult.length == 2) {
            orgin = cuppaResult.split(" \\(")[0];
            prediction = cuppaResult.split("\\(")[1];
            prediction = prediction.substring(0, prediction.length() - 1);
        } else if (lengthCuppaResult.length == 1) {
            orgin = cuppaResult.split("\\(")[0];
        } else {
            LOGGER.warn("Cuppa result is Empty");
        }

        return ImmutableMolecularTissueOrginData.builder()
                .conclusion(cuppaResult)
                .predictedOrigin(orgin)
                .predictionLikelihood(prediction)
                .build();
    }
}