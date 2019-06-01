package com.hartwig.hmftools.common.variant.structural.linx;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class LinxClusterFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";
    private static final String FILE_EXTENSION = ".linx_cluster.csv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<LinxSvData> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<LinxSvData> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<LinxSvData> svDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        svDataList.stream().map(LinxClusterFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<LinxSvData> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().map(LinxClusterFile::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER)
                .add("svId")
                .add("clusterId")
                .add("clusterReason")
                .add("fragileSiteStart")
                .add("fragileSiteEnd")
                .add("isFoldback")
                .add("lineTypeStart")
                .add("lineTypeEnd")
                .add("ploidyMin")
                .add("ploidyMax")
                .add("geneStart")
                .add("geneEnd")
                .add("replicationTimingStart")
                .add("replicationTimingEnd")
                .add("localTopologyIdStart")
                .add("localTopologyIdEnd")
                .add("localTopologyStart")
                .add("localTopologyEnd")
                .add("localTICountStart")
                .add("localTICountEnd")
                .add("nearestTypeStart")
                .add("nearestTypeEnd")
                .add("nearestLengthStart")
                .add("nearestLengthEnd")
                .add("dbLengthStart")
                .add("dbLengthEnd")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final LinxSvData svData) {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(svData.svId()))
                .add(String.valueOf(svData.clusterId()))
                .add(String.valueOf(svData.clusterReason()))
                .add(String.valueOf(svData.fragileSiteStart()))
                .add(String.valueOf(svData.fragileSiteEnd()))
                .add(String.valueOf(svData.isFoldback()))
                .add(String.valueOf(svData.lineTypeStart()))
                .add(String.valueOf(svData.lineTypeEnd()))
                .add(FORMAT.format(svData.ploidyMin()))
                .add(FORMAT.format(svData.ploidyMax()))
                .add(String.valueOf(svData.geneStart()))
                .add(String.valueOf(svData.geneEnd()))
                .add(FORMAT.format(svData.replicationTimingStart()))
                .add(FORMAT.format(svData.replicationTimingEnd()))
                .add(String.valueOf(svData.localTopologyIdStart()))
                .add(String.valueOf(svData.localTopologyIdEnd()))
                .add(String.valueOf(svData.localTopologyStart()))
                .add(String.valueOf(svData.localTopologyEnd()))
                .add(String.valueOf(svData.localTICountStart()))
                .add(String.valueOf(svData.localTICountEnd()))
                .add(String.valueOf(svData.nearestTypeStart()))
                .add(String.valueOf(svData.nearestTypeEnd()))
                .add(String.valueOf(svData.nearestLengthStart()))
                .add(String.valueOf(svData.nearestLengthEnd()))
                .add(String.valueOf(svData.dbLengthStart()))
                .add(String.valueOf(svData.dbLengthEnd()))
                .toString();
    }

    @NotNull
    private static LinxSvData fromString(@NotNull final String clusterData)
    {
        String[] values = clusterData.split(DELIMITER);

        int index = 0;

        return ImmutableLinxSvData.builder()
                .svId(Integer.valueOf(values[index++]))
                .clusterId(Integer.valueOf(values[index++]))
                .clusterReason(values[index++])
                .fragileSiteStart(Boolean.valueOf(values[index++]))
                .fragileSiteEnd(Boolean.valueOf(values[index++]))
                .isFoldback(Boolean.valueOf(values[index++]))
                .lineTypeStart(values[index++])
                .lineTypeEnd(values[index++])
                .ploidyMin(Double.valueOf(values[index++]))
                .ploidyMax(Double.valueOf(values[index++]))
                .geneStart(values[index++])
                .geneEnd(values[index++])
                .replicationTimingStart(Double.valueOf(values[index++]))
                .replicationTimingEnd(Double.valueOf(values[index++]))
                .localTopologyIdStart(Integer.valueOf(values[index++]))
                .localTopologyIdEnd(Integer.valueOf(values[index++]))
                .localTopologyStart(values[index++])
                .localTopologyEnd(values[index++])
                .localTICountStart(Integer.valueOf(values[index++]))
                .localTICountEnd(Integer.valueOf(values[index++]))
                .nearestTypeStart(values[index++])
                .nearestTypeEnd(values[index++])
                .nearestLengthStart(Long.valueOf(values[index++]))
                .nearestLengthEnd(Long.valueOf(values[index++]))
                .dbLengthStart(Long.valueOf(values[index++]))
                .dbLengthEnd(Long.valueOf(values[index++]))
                .build();
    }
}
