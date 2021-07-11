package com.hartwig.hmftools.common.sv.linx;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

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
public abstract class LinxCluster
{
    public abstract int clusterId();
    public abstract String category();
    public abstract boolean synthetic();
    public abstract String resolvedType();
    public abstract int clusterCount();
    public abstract String clusterDesc();

    public static final String DELIMITER = "\t";

    private static final String FILE_EXTENSION = ".linx.clusters.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<LinxCluster> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<LinxCluster> clusters) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(clusters));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<LinxCluster> clusters)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        clusters.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<LinxCluster> fromLines(@NotNull List<String> lines)
    {
        final String header = lines.get(0);
        lines.remove(0);

        if(header.contains("subClonal"))
        {
            final Map<String,Integer> fieldIndexMap = createFieldsIndexMap(header,DELIMITER);
            return lines.stream().map(x -> fromString_v1_10(x, fieldIndexMap)).collect(toList());
        }
        else
        {
            return lines.stream().map(LinxCluster::fromString).collect(toList());
        }
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("clusterId")
                .add("category")
                .add("synthetic")
                .add("resolvedType")
                .add("clusterCount")
                .add("clusterDesc")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final LinxCluster cluster)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(cluster.clusterId()))
                .add(String.valueOf(cluster.category()))
                .add(String.valueOf(cluster.synthetic()))
                .add(String.valueOf(cluster.resolvedType()))
                .add(String.valueOf(cluster.clusterCount()))
                .add(String.valueOf(cluster.clusterDesc()))
                .toString();
    }

    @NotNull
    private static LinxCluster fromString(@NotNull final String clusterData)
    {
        String[] values = clusterData.split(DELIMITER);

        int index = 0;

        return ImmutableLinxCluster.builder()
                .clusterId(Integer.parseInt(values[index++]))
                .category(values[index++])
                .synthetic(Boolean.parseBoolean(values[index++]))
                .resolvedType(values[index++])
                .clusterCount(Integer.parseInt(values[index++]))
                .clusterDesc(values[index++])
                .build();
    }

    @NotNull
    private static LinxCluster fromString_v1_10(@NotNull final String clusterData, final Map<String,Integer> fieldIndexMap)
    {
        String[] values = clusterData.split(DELIMITER);

        return ImmutableLinxCluster.builder()
                .clusterId(Integer.parseInt(values[fieldIndexMap.get("clusterId")]))
                .category(values[fieldIndexMap.get("resolvedType")])
                .synthetic(Boolean.parseBoolean(values[fieldIndexMap.get("synthetic")]))
                .resolvedType(values[fieldIndexMap.get("subType")])
                .clusterCount(Integer.parseInt(values[fieldIndexMap.get("clusterCount")]))
                .clusterDesc(values[fieldIndexMap.get("clusterDesc")])
                .build();
    }

}
