package com.hartwig.hmftools.common.sv.linx;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.sv.linx.LinxCluster.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
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
    private static final String CATALOG_EXTENSION = ".linx.driver.catalog.tsv";
    private static final String GERMLINE_CATALOG_EXTENSION = ".linx.germline.driver.catalog.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static String generateCatalogFilename(@NotNull final String basePath, @NotNull final String sample, boolean isSomatic)
    {
        if(isSomatic)
            return basePath + File.separator + sample + CATALOG_EXTENSION;
        else
            return basePath + File.separator + sample + GERMLINE_CATALOG_EXTENSION;
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

    private static List<String> toLines(final List<LinxDriver> drivers)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        drivers.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<LinxDriver> fromLines(final List<String> lines)
    {
        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);
        lines.remove(0);

        List<LinxDriver> drivers = Lists.newArrayList();

        for(int i = 0; i < lines.size(); ++i)
        {
            String[] values = lines.get(i).split(DELIMITER);

            drivers.add(ImmutableLinxDriver.builder()
                    .clusterId(Integer.parseInt(values[fieldsIndexMap.get("clusterId")]))
                    .gene(values[fieldsIndexMap.get("gene")])
                    .eventType(values[fieldsIndexMap.get("eventType")])
                    .build());
        }

        return drivers;
    }

    private static String header()
    {
        return new StringJoiner(LinxCluster.DELIMITER)
                .add("clusterId")
                .add("gene")
                .add("eventType")
                .toString();
    }

    private static String toString(@NotNull final LinxDriver driver)
    {
        return new StringJoiner(LinxCluster.DELIMITER)
                .add(String.valueOf(driver.clusterId()))
                .add(String.valueOf(driver.gene()))
                .add(String.valueOf(driver.eventType()))
                .toString();
    }
}
