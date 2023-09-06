package com.hartwig.hmftools.common.circos;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

import org.jetbrains.annotations.NotNull;

public final class CircosINDELWriter
{
    public static void writePositions(final String filePath, Collection<VariantContextDecorator> values)
            throws IOException
    {
        writeCircosFile(filePath, values, CircosINDELWriter::transformPosition);
    }

    private static <T> void writeCircosFile(final String filePath, Collection<T> values, Function<T, String> toStringFunction) throws IOException
    {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        values.stream().map(toStringFunction).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    private static String header()
    {
        return "#chromosome\tstart\tend\tvalue";
    }

    private static String transformPosition(final VariantContextDecorator position)
    {
        return CircosFileWriter.transformPosition(position, CircosINDELWriter::color);
    }

    private static String color(final VariantContextDecorator variant)
    {
        if(variant.alt().length() > variant.ref().length())
        {
            return "vdyellow";
        }

        return "red";
    }
}
