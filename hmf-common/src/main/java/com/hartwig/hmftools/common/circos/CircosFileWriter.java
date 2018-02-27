package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class CircosFileWriter {

    public static <T extends GenomeRegion> void writeRegions(@NotNull final String filePath,
            @NotNull Collection<T> values, @NotNull ToDoubleFunction<T> valueExtractor) throws IOException {
        Function<T, String> toString = t -> transformRegion(valueExtractor, t);
        writeCircosFile(filePath, values, toString);
    }

    public static <T extends GenomePosition> void writePositions(@NotNull final String filePath,
            @NotNull Collection<T> values, @NotNull ToDoubleFunction<T> valueExtractor) throws IOException {
        Function<T, String> toString = t -> transformPosition(valueExtractor, t);
        writeCircosFile(filePath, values, toString);
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

    private static <T extends GenomePosition> String transformPosition(ToDoubleFunction<T> valueExtractor,
            T position) {
        return new StringBuilder().append("hs")
                .append(position.chromosome())
                .append('\t')
                .append(position.position())
                .append('\t')
                .append(position.position())
                .append('\t')
                .append(valueExtractor.applyAsDouble(position))
                .toString();
    }

    private static <T extends GenomeRegion> String transformRegion(ToDoubleFunction<T> valueExtractor, T region) {
        return new StringBuilder().append("hs")
                .append(region.chromosome())
                .append('\t')
                .append(region.start())
                .append('\t')
                .append(region.end())
                .append('\t')
                .append(valueExtractor.applyAsDouble(region))
                .toString();
    }

}
