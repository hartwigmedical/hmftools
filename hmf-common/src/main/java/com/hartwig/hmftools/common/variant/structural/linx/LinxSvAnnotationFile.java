package com.hartwig.hmftools.common.variant.structural.linx;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.variant.structural.linx.LinxClusterFile.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class LinxSvAnnotationFile
{
    private static final String FILE_EXTENSION = ".linx.svs.tsv";

    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.0000");

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<LinxSvAnnotation> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<LinxSvAnnotation> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<LinxSvAnnotation> svDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        svDataList.stream().map(LinxSvAnnotationFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<LinxSvAnnotation> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("svId")).map(LinxSvAnnotationFile::fromString).collect(toList());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("svId")
                .add("clusterId")
                .add("clusterReason")
                .add("fragileSiteStart")
                .add("fragileSiteEnd")
                .add("isFoldback")
                .add("lineTypeStart")
                .add("lineTypeEnd")
                .add("junctionCopyNumberMin")
                .add("junctionCopyNumberMax")
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
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final LinxSvAnnotation svData) {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(svData.svId()))
                .add(String.valueOf(svData.clusterId()))
                .add(String.valueOf(svData.clusterReason()))
                .add(String.valueOf(svData.fragileSiteStart()))
                .add(String.valueOf(svData.fragileSiteEnd()))
                .add(String.valueOf(svData.isFoldback()))
                .add(String.valueOf(svData.lineTypeStart()))
                .add(String.valueOf(svData.lineTypeEnd()))
                .add(DECIMAL_FORMAT.format(svData.junctionCopyNumberMin()))
                .add(DECIMAL_FORMAT.format(svData.junctionCopyNumberMax()))
                .add(String.valueOf(svData.geneStart()))
                .add(String.valueOf(svData.geneEnd()))
                .add(DECIMAL_FORMAT.format(svData.replicationTimingStart()))
                .add(DECIMAL_FORMAT.format(svData.replicationTimingEnd()))
                .add(String.valueOf(svData.localTopologyIdStart()))
                .add(String.valueOf(svData.localTopologyIdEnd()))
                .add(String.valueOf(svData.localTopologyStart()))
                .add(String.valueOf(svData.localTopologyEnd()))
                .add(String.valueOf(svData.localTICountStart()))
                .add(String.valueOf(svData.localTICountEnd()))
                .toString();
    }

    @NotNull
    private static LinxSvAnnotation fromString(@NotNull final String svData)
    {
        String[] values = svData.split(DELIMITER);

        int index = 0;

        return ImmutableLinxSvAnnotation.builder()
                .svId(Integer.parseInt(values[index++]))
                .clusterId(Integer.parseInt(values[index++]))
                .clusterReason(values[index++])
                .fragileSiteStart(Boolean.parseBoolean(values[index++]))
                .fragileSiteEnd(Boolean.parseBoolean(values[index++]))
                .isFoldback(Boolean.parseBoolean(values[index++]))
                .lineTypeStart(values[index++])
                .lineTypeEnd(values[index++])
                .junctionCopyNumberMin(Double.parseDouble(values[index++]))
                .junctionCopyNumberMax(Double.parseDouble(values[index++]))
                .geneStart(values[index++])
                .geneEnd(values[index++])
                .replicationTimingStart(Double.parseDouble(values[index++]))
                .replicationTimingEnd(Double.parseDouble(values[index++]))
                .localTopologyIdStart(Integer.parseInt(values[index++]))
                .localTopologyIdEnd(Integer.parseInt(values[index++]))
                .localTopologyStart(values[index++])
                .localTopologyEnd(values[index++])
                .localTICountStart(Integer.parseInt(values[index++]))
                .localTICountEnd(Integer.parseInt(values[index++]))
                .build();
    }
}
