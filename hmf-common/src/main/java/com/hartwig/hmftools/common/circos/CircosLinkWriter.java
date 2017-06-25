package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.StringJoiner;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

public class CircosLinkWriter {

    public static void writeVariants(@NotNull final String filePath, @NotNull Collection<StructuralVariant> values)
            throws IOException {
        writeCircosFile(filePath, values, CircosLinkWriter::toString);
    }

    private static <T> void writeCircosFile(@NotNull final String filePath, @NotNull Collection<T> values,
            Function<T, String> toStringFunction) throws IOException {
        final Collection<String> lines = Lists.newArrayList();
        values.stream().map(toStringFunction).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    private static String toString(StructuralVariant variant) {
        return new StringJoiner("\t").add("hs" + variant.startChromosome())
                .add(String.valueOf(variant.startPosition()))
                .add(String.valueOf(variant.startPosition()))
                .add("hs" + variant.endChromosome())
                .add(String.valueOf(variant.endPosition()))
                .add(String.valueOf(variant.endPosition()))
                .add("color=" + color(variant))
                .toString();
    }

    private static String color(StructuralVariant variant) {
        switch (variant.type()) {
            case DUP:
                return "green";
            case DEL:
                return "red";
            case BND:
                return "blue";
            case INS:
                return "vdyellow";
            case INV:
                return "black";
            default:
                return "purple";
        }
    }
}
