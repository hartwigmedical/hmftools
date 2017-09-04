package com.hartwig.hmftools.common.ecrf.doid;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class TumorLocationDoidMapping {

    private static final Logger LOGGER = LogManager.getLogger(TumorLocationDoidMapping.class);

    private static final int TUMOR_TYPE_COLUMN = 0;
    private static final int DOID_COLUMN = 1;
    private static final int FIELD_COUNT = 2;

    private static final String FIELD_SEPARATOR = ",";

    private static final Pattern DOID_PATTERN = Pattern.compile("DOID:[0-9]+");

    private TumorLocationDoidMapping() {
    }

    @NotNull
    public static Map<String, Set<String>> readMappingFromCSV(@NotNull final String pathToCsv) throws IOException, EmptyFileException {
        final Map<String, Set<String>> recipientsPerCentra = Maps.newHashMap();
        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT);
            if (parts.length == FIELD_COUNT) {
                recipientsPerCentra.put(parts[TUMOR_TYPE_COLUMN], extractDoids(parts[DOID_COLUMN]));
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line: " + line);
            }
        }
        return recipientsPerCentra;
    }

    @NotNull
    private static Set<String> extractDoids(@NotNull final String doidString) {
        final Matcher matcher = DOID_PATTERN.matcher(doidString);
        final Set<String> result = Sets.newHashSet();
        while (matcher.find()) {
            result.add(matcher.group(0));
        }
        return result;
    }
}
