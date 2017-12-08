package com.hartwig.hmftools.common.ecrf.doid;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TumorLocationDoidMapping {

    private static final Logger LOGGER = LogManager.getLogger(TumorLocationDoidMapping.class);

    private static final int TUMOR_TYPE_COLUMN = 0;
    private static final int DOID_COLUMN = 1;
    private static final int FIELD_COUNT = 2;

    private static final String FIELD_SEPARATOR = ",";

    private static final Pattern DOID_PATTERN = Pattern.compile("DOID:[0-9]+");

    @NotNull
    protected abstract Map<String, Set<String>> doidsPerTumorType();

    @NotNull
    public Set<String> doidsForTumorType(@NotNull final String tumorType) {
        final Set<String> result = doidsPerTumorType().get(tumorType);
        return result == null ? Sets.newHashSet() : result;
    }

    @NotNull
    static TumorLocationDoidMapping fromCSV(@NotNull final String pathToCsv) throws IOException, EmptyFileException {
        return fromLines(FileReader.build().readLines(new File(pathToCsv).toPath()));
    }

    @NotNull
    public static TumorLocationDoidMapping fromResource(@NotNull final String name) {
        final InputStream inputStream = TumorLocationDoidMapping.class.getResourceAsStream(name);
        return fromLines(new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static TumorLocationDoidMapping fromLines(@NotNull final List<String> lines) {
        final Map<String, Set<String>> doidsPerTumorType = Maps.newHashMap();
        for (String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT);
            if (parts.length == FIELD_COUNT) {
                doidsPerTumorType.put(parts[TUMOR_TYPE_COLUMN], extractDoids(parts[DOID_COLUMN]));
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line: " + line);
            }
        }
        return ImmutableTumorLocationDoidMapping.of(doidsPerTumorType);
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
