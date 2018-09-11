package com.hartwig.hmftools.actionability.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class ActionabilityAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityAnalyzer.class);
    static final String DELIMITER = "\t";

    @NotNull
    public static ActionabilityAnalyzer loadFromFile(String file) throws IOException {
        final String line =  Files.readAllLines(new File(file).toPath()).get(1);

        return null;
    }



}
