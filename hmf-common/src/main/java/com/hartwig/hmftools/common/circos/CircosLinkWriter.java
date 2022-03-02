package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.jetbrains.annotations.NotNull;

public final class CircosLinkWriter
{
    private CircosLinkWriter()
    {
    }

    public static void writeVariants(@NotNull final String filePath, @NotNull Collection<StructuralVariant> values) throws IOException
    {
        final Collection<String> lines = values.stream()
                .filter(x -> x.end() != null)
                .map(CircosLinkWriter::toString)
                .collect(Collectors.toList());
        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    private static String toString(@NotNull final StructuralVariant variant)
    {
        String startChromosome = variant.chromosome(true);
        String endChromosome = variant.chromosome(false);

        // SGLs are filtered out at this stage
        assert startChromosome != null && endChromosome != null;

        return new StringJoiner("\t").add(CircosFileWriter.circosContig(startChromosome))
                .add(String.valueOf(variant.position(true)))
                .add(String.valueOf(variant.position(true)))
                .add(CircosFileWriter.circosContig(endChromosome))
                .add(String.valueOf(variant.position(false)))
                .add(String.valueOf(variant.position(false)))
                .add("color=" + color(variant))
                .toString();
    }

    @NotNull
    private static String color(@NotNull final StructuralVariant variant)
    {
        switch(variant.type())
        {
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
