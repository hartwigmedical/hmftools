package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.StringJoiner;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class CircosINDELWriter {

    public static void writePositions(@NotNull final String filePath, @NotNull Collection<PurityAdjustedSomaticVariant> values)
            throws IOException {
        writeCircosFile(filePath, values, CircosINDELWriter::transformPosition);
    }

    private static <T> void writeCircosFile(@NotNull final String filePath, @NotNull Collection<T> values,
            Function<T, String> toStringFunction) throws IOException {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        values.stream().map(toStringFunction).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    private static String header() {
        return "#chromosome\tstart\tend\tvalue";
    }

    private static String transformPosition(PurityAdjustedSomaticVariant position) {
        return new StringJoiner("\t").add("hs" + position.chromosome())
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.adjustedVAF()))
                .add("color=" + color(position))
                .toString();
    }

    private static String color(PurityAdjustedSomaticVariant variant) {

        if (variant.alt().length() > variant.ref().length()) {
            return "vdyellow";
        }

        return "red";
    }

}
