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
    private CircosINDELWriter()
    {
    }

    public static void writePositions(@NotNull final String filePath, @NotNull Collection<VariantContextDecorator> values)
            throws IOException
    {
        writeCircosFile(filePath, values, CircosINDELWriter::transformPosition);
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
        return CircosFileWriter.transformPosition(position, CircosINDELWriter::color);
    }

    @NotNull
    private static String color(@NotNull final VariantContextDecorator variant)
    {
        if(variant.alt().length() > variant.ref().length())
        {
            return "vdyellow";
        }

        return "red";
    }
}
