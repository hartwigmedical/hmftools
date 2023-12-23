package com.hartwig.hmftools.esvee;

import com.google.errorprone.annotations.MustBeClosed;
import com.hartwig.hmftools.common.sv.Direction;
import com.hartwig.hmftools.esvee.util.StringUtils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Objects;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Stream;

public enum JunctionReader {
    ;

    @MustBeClosed
    public static Stream<Junction> readJunctionFile(final File junctionFile) {
        try {
            //noinspection resource
            return Files.lines(junctionFile.toPath()).map(new Function<String, Junction>() {
                boolean hadHeader = false;
                String delimiter = ",";

                @Override
                public Junction apply(final String line) {
                    if (!hadHeader) {
                        hadHeader = true;
                        delimiter = guessDelimiter(line);
                        return null;
                    }

                    return parseJunction(line, delimiter);
                }
            }).filter(Objects::nonNull);
        } catch (final IOException exception) {
            throw new RuntimeException(exception);
        }
    }

    private static String guessDelimiter(final String firstLine) {
        if (!firstLine.toLowerCase().startsWith("chromosome")) {
            throw new IllegalArgumentException("Line does not appear to be the start of a junction file: " + firstLine);
        }

        final String delimiter = firstLine.substring("chromosome".length(), "chromosome".length() + 1);
        if (firstLine.split(Pattern.quote(delimiter)).length != 15) {
            throw new IllegalArgumentException("Expected to find 15 fields in header (using delim '" + delimiter + "'), found: " + firstLine);
        }

        return delimiter;
    }

    private static Junction parseJunction(final String line, final String delimiter) {
        final String[] parts = line.split(Pattern.quote(delimiter));
        if (parts.length != 15) {
            throw new IllegalArgumentException("Line does not have the correct number of fields");
        }

        return ImmutableJunction.builder()
                .chromosome(parts[0])
                .position(Integer.parseInt(parts[1]))
                .orientation(Integer.parseInt(parts[2]) > 0 ? Direction.FORWARDS : Direction.REVERSE)
                .junctionFragments(Integer.parseInt(parts[3]))
                .supportFragments(Integer.parseInt(parts[4]))
                .discordantFragments(Integer.parseInt(parts[5]))
                .lowMapQualityFragments(Integer.parseInt(parts[6]))
                .maxMapQuality(Integer.parseInt(parts[7]))
                .maxSoftClipLength(Integer.parseInt(parts[8]))
                .baseDepth(Integer.parseInt(parts[9]))
                .hasPolyAT(StringUtils.parseBoolean(parts[10]))
                .isIndel(StringUtils.parseBoolean(parts[11]))
                .isHotspot(StringUtils.parseBoolean(parts[12]))
                .softClipBases(parts[13])
                .initialReadId(parts[14])
                .build();
    }
}
