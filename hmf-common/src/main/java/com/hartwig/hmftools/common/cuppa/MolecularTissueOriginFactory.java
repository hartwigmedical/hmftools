package com.hartwig.hmftools.common.cuppa;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class MolecularTissueOriginFactory {

    private MolecularTissueOriginFactory(){
    }

    private static final Logger LOGGER = LogManager.getLogger(MolecularTissueOriginFactory.class);


    @NotNull
    public static String readMolecularTissueOriginResult(@NotNull String molecularTissueOriginTxt) throws IOException {
        return fromLines(Files.readAllLines(new File(molecularTissueOriginTxt).toPath()));
    }

    @NotNull
    static String fromLines(@NotNull final List<String> lines) {
        if (lines.size() == 1) {
            return lines.get(0).split(" - ")[1];
        } else {
            LOGGER.warn("No molecular tissue origin is known!");
            return Strings.EMPTY;
        }

    }
}
