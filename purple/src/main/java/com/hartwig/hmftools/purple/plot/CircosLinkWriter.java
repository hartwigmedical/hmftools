package com.hartwig.hmftools.purple.plot;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CircosLinkWriter
{
    public static void writeVariants(final String filePath, Collection<StructuralVariant> values) throws IOException
    {
        final Collection<String> lines = values.stream()
                .filter(x -> x.end() != null)
                .map(CircosLinkWriter::toString)
                .collect(Collectors.toList());
        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    private static String toString(final StructuralVariant variant)
    {
        String startChromosome = variant.chromosome(true);
        String endChromosome = variant.chromosome(false);

        // SGLs are filtered out at this stage
        assert startChromosome != null && endChromosome != null;

        String bezierRadiusVal = bezierRadius(variant);

        StringJoiner stringJoiner = new StringJoiner("\t").add(CircosFileWriter.circosContig(startChromosome))
                .add(String.valueOf(variant.position(true)))
                .add(String.valueOf(variant.position(true)))
                .add(CircosFileWriter.circosContig(endChromosome))
                .add(String.valueOf(variant.position(false)))
                .add(String.valueOf(variant.position(false)));

        if(bezierRadiusVal == null)
        {
            stringJoiner.add("color=" + color(variant));
        }
        else
        {
            stringJoiner.add("color=" + color(variant) + ",bezier_radius=" + bezierRadiusVal);
        }

        return stringJoiner.toString();
    }

    @NotNull
    private static String color(final StructuralVariant variant)
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

    // calculate the bezier radius of the link.
    // it is proportional to the log10 size of the SV if it is within same chromosome
    @Nullable
    static String bezierRadius(final StructuralVariant variant)
    {
        String startChromosome = variant.chromosome(true);
        String endChromosome = variant.chromosome(false);

        // SGLs are filtered out at this stage
        if (startChromosome == null || !startChromosome.equals(endChromosome))
        {
            return null;
        }

        int length;

        Integer start = variant.position(true);
        Integer end = variant.position(false);

        if(start != null && end != null)
        {
            // end is inclusive
            length = Math.max(end - start + 1, 1);
        }
        else
        {
            return null;
        }

        // between 0.1 to 0.35
        // which is scaled by log scale
        // we set 10m to be max, which is 1e7
        // the bez radius would be 0.1 when length is 10m, and 0.35 when length is 1

        double maxBezRadius = 0.35; // this value is the same as in the circos.template link section (radius)
        double minBezRadius = 0.1; // this value is the same as in the circos.template link section (bezier radius)
        double log10Max = 7.0;
        double r = Math.min(Math.log10(length), log10Max);

        r = maxBezRadius - r / log10Max * (maxBezRadius - minBezRadius);

        return String.format("%.3fr", r);
    }
}
