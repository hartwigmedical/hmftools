package com.hartwig.hmftools.common.flagstat;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FlagstatFile {

    private FlagstatFile() {
    }

    @NotNull
    public static Flagstat read(@NotNull String flagstatsPath) throws IOException {
        List<String> lines = Files.readAllLines(new File(flagstatsPath).toPath());
        String total = valueBySubstring(lines, "total");
        String mapped = valueBySubstring(lines, "mapped (");
        String duplicates = valueBySubstring(lines, "duplicates");
        if (total == null || mapped == null || duplicates == null) {
            throw new IOException("Unable to parse flagstat file correctly");
        }

        return ImmutableFlagstat.builder()
                .mappedProportion(divideTwoStrings(mapped, total))
                .duplicateProportion(divideTwoStrings(duplicates, total))
                .build();
    }

    private static double divideTwoStrings(@NotNull String string1, @NotNull String string2) {
        long value1 = Long.parseLong(string1);
        long value2 = Long.parseLong(string2);
        return roundToSixDecimals((double) value1 / value2);
    }

    private static double roundToSixDecimals(double value) {
        BigDecimal bd = BigDecimal.valueOf(value);
        bd = bd.setScale(6, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

    @Nullable
    private static String valueBySubstring(@NotNull List<String> lines, @NotNull String subString) {
        List<String> matchLines = Lists.newArrayList();
        for (String line : lines) {
            if (line.contains(subString)) {
                matchLines.add(line);
            }
        }
        if (matchLines.size() == 1) {
            return matchLines.get(0).split(" ")[0];
        }
        return null;
    }
}
