package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

import org.jetbrains.annotations.NotNull;

public final class CircosSNPWriter
{
    private static final String PINK = "233,187,184";
    private static final String BLUE = "20,176,239";
    private static final String BLACK = "6,8,9";
    private static final String RED = "224,7,20";
    private static final String GREY = "191,190,191";
    private static final String GREEN = "144,202,75";

    private CircosSNPWriter()
    {
    }

    public static void writePositions(@NotNull final String filePath, @NotNull Collection<VariantContextDecorator> values)
            throws IOException
    {
        writeCircosFile(filePath, values, CircosSNPWriter::transformPosition);
    }

    private static <T> void writeCircosFile(@NotNull final String filePath, @NotNull Collection<T> values,
            @NotNull Function<T, String> toStringFunction) throws IOException
    {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        values.stream().map(toStringFunction).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    private static String header()
    {
        return "#chromosome\tstart\tend\tvalue";
    }

    @NotNull
    private static String transformPosition(@NotNull final VariantContextDecorator position)
    {
        return CircosFileWriter.transformPosition(position, CircosSNPWriter::color);
    }

    @NotNull
    private static String color(@NotNull final VariantContextDecorator variant)
    {
        if(signature("C", "A", variant))
        {
            return BLUE;
        }
        if(signature("C", "G", variant))
        {
            return BLACK;
        }
        if(signature("C", "T", variant))
        {
            return RED;
        }
        if(signature("T", "A", variant))
        {
            return GREY;
        }
        if(signature("T", "C", variant))
        {
            return GREEN;
        }
        if(signature("T", "G", variant))
        {
            return PINK;
        }

        return "purple";
    }

    private static boolean signature(@NotNull final String ref, @NotNull final String alt,
            @NotNull final VariantContextDecorator variant)
    {
        return (variant.ref().equals(ref) && variant.alt().equals(alt)) || (variant.ref().equals(inverse(ref)) && variant.alt()
                .equals(inverse(alt)));
    }

    @NotNull
    private static String inverse(@NotNull final String base)
    {
        if(base.equals("G"))
        {
            return "C";
        }
        if(base.equals("C"))
        {
            return "G";
        }
        if(base.equals("A"))
        {
            return "T";
        }
        return "A";
    }
}
