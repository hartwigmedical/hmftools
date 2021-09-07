package com.hartwig.hmftools.common.drivercatalog.dnds;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class DndsMutationalLoadFile
{
    private static final String DELIMITER = "\t";

    private DndsMutationalLoadFile()
    {
    }

    @NotNull
    public static List<DndsMutationalLoad> read(@NotNull final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull final List<DndsMutationalLoad> variants) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(true, variants));
    }

    public static void writeHeader(@NotNull final String filename) throws IOException
    {
        Files.write(new File(filename).toPath(), Collections.singletonList(header()));
    }

    public static void append(@NotNull final String filename, @NotNull final List<DndsMutationalLoad> load) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(false, load), StandardOpenOption.APPEND);
    }

    @NotNull
    private static List<DndsMutationalLoad> fromLines(@NotNull final List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("sample")).map(DndsMutationalLoadFile::fromString).collect(Collectors.toList());
    }

    @NotNull
    private static List<String> toLines(boolean header, @NotNull final List<DndsMutationalLoad> variants)
    {
        final List<String> lines = Lists.newArrayList();
        if(header)
        {
            lines.add(header());
        }
        variants.stream().map(DndsMutationalLoadFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("sample")
                .add("snvBiallelic")
                .add("snvNonBiallelic")
                .add("indelBiallelic")
                .add("indelNonBiallelic")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final DndsMutationalLoad variant)
    {
        return new StringJoiner(DELIMITER).add(variant.sampleId())
                .add(String.valueOf(variant.snvBiallelic()))
                .add(String.valueOf(variant.snvNonBiallelic()))
                .add(String.valueOf(variant.indelBiallelic()))
                .add(String.valueOf(variant.indelNonBiallelic()))
                .toString();
    }

    @NotNull
    private static DndsMutationalLoad fromString(@NotNull final String line)
    {
        String[] values = line.split(DELIMITER);
        return ImmutableDndsMutationalLoad.builder()
                .sampleId(values[0])
                .snvBiallelic(Integer.parseInt(values[1]))
                .snvNonBiallelic(Integer.parseInt(values[2]))
                .indelBiallelic(Integer.parseInt(values[3]))
                .indelNonBiallelic(Integer.parseInt(values[4]))
                .build();
    }
}
