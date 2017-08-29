package com.hartwig.hmftools.common.center;

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

public final class Center {

    private static final Logger LOGGER = LogManager.getLogger(Center.class);

    private static final int ID_COLUMN = 0;
    private static final int CPCT_RECIPIENTS_COLUMN = 4;
    private static final int DRUP_RECIPIENTS_COLUMN = 7;
    private static final int FIELD_COUNT = 13;
    private static final int ADDRESS_NAME_FIELD = 10;
    private static final int ADDRESS_ZIP_FIELD = 11;
    private static final int ADDRESS_CITY_FIELD = 12;

    private static final String FIELD_SEPARATOR = ",";

    private Center() {
    }

    @NotNull
    public static CenterModel readFromCSV(@NotNull final String pathToCsv) throws IOException, EmptyFileException {
        final Map<String, CenterData> centerPerId = Maps.newHashMap();
        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (final String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT);
            if (parts.length == FIELD_COUNT) {
                centerPerId.put(parts[ID_COLUMN],
                        ImmutableCenterData.of(parts[CPCT_RECIPIENTS_COLUMN], parts[DRUP_RECIPIENTS_COLUMN], parts[ADDRESS_NAME_FIELD],
                                parts[ADDRESS_ZIP_FIELD], parts[ADDRESS_CITY_FIELD]));
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in center csv: " + line);
            }
        }
        return ImmutableCenterModel.of(centerPerId);
    }
}
