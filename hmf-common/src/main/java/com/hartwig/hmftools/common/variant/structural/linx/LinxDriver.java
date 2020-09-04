package com.hartwig.hmftools.common.variant.structural.linx;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.variant.structural.linx.LinxCluster.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
public abstract class LinxDriver
{
    public abstract int clusterId();
    public abstract String gene();
    public abstract String eventType();

    private static final String FILE_EXTENSION = ".linx.drivers.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<LinxDriver> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<LinxDriver> clusters) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(clusters));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<LinxDriver> drivers)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        drivers.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<LinxDriver> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("clusterId")).map(LinxDriver::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER)
                .add("clusterId")
                .add("gene")
                .add("eventType")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final LinxDriver driver)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(driver.clusterId()))
                .add(String.valueOf(driver.gene()))
                .add(String.valueOf(driver.eventType()))
                .toString();
    }

    @NotNull
    private static LinxDriver fromString(@NotNull final String clusterData)
    {
        String[] values = clusterData.split(DELIMITER);

        int index = 0;

        return ImmutableLinxDriver.builder()
                .clusterId(Integer.parseInt(values[index++]))
                .gene(values[index++])
                .eventType(values[index++])
                .build();
    }


}
