package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

public class CircosLinkWriter {

    public static void writeVariants(@NotNull final String filePath, @NotNull Collection<StructuralVariant> values) throws IOException {
        final Collection<String> lines = values.stream()
                .filter(x -> x.end() != null)
                .map(CircosLinkWriter::toString)
                .collect(Collectors.toList());
        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    private static String toString(@NotNull final StructuralVariant variant) {
        return new StringJoiner("\t").add(CircosFileWriter.circosContig(variant.chromosome(true)))
                .add(String.valueOf(variant.position(true)))
                .add(String.valueOf(variant.position(true)))
                .add(CircosFileWriter.circosContig(variant.chromosome(false)))
                .add(String.valueOf(variant.position(false)))
                .add(String.valueOf(variant.position(false)))
                .add("color=" + color(variant))
                .toString();
    }

    @NotNull
    private static String color(@NotNull final StructuralVariant variant) {
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
