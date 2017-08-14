package com.hartwig.hmftools.common.centra;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class Centra {

    private static final Logger LOGGER = LogManager.getLogger(Centra.class);

    private static final int ID_COLUMN = 0;
    private static final int CPCT_RECIPIENTS_COLUMN = 4;
    private static final int DRUP_RECIPIENTS_COLUMN = 7;
    private static final int FIELD_COUNT = 11;

    private static final String FIELD_SEPARATOR = ",";

    private Centra() {
    }

    @NotNull
    public static Map<String, String> readCPCTRecipientsFromCSV(@NotNull final String pathToCsv) throws IOException, EmptyFileException {
        final Map<String, String> recipientsPerCentra = Maps.newHashMap();
        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT);
            if (parts.length == FIELD_COUNT) {
                recipientsPerCentra.put(parts[ID_COLUMN], parts[CPCT_RECIPIENTS_COLUMN]);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in form status csv: " + line);
            }
        }
        return recipientsPerCentra;
    }

    @NotNull
    public static Map<String, String> readDRUPRecipientsFromCSV(@NotNull final String pathToCsv) throws IOException, EmptyFileException {
        final Map<String, String> recipientsPerCentra = Maps.newHashMap();
        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT);
            if (parts.length == FIELD_COUNT) {
                recipientsPerCentra.put(parts[ID_COLUMN], parts[DRUP_RECIPIENTS_COLUMN]);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in form status csv: " + line);
            }
        }
        return recipientsPerCentra;
    }
}
