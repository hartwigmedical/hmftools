package com.hartwig.hmftools.common.variant.structural.linx;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.variant.structural.linx.LinxClusterFile.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class LinxLinkFile
{
    private static final String FILE_EXTENSION = ".linx.links.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<LinxLink> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<LinxLink> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<LinxLink> svDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        svDataList.stream().map(LinxLinkFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<LinxLink> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("clusterId")).map(LinxLinkFile::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER)
                .add("clusterId")
                .add("chainId")
                .add("chainIndex")
                .add("chainCount")
                .add("lowerBreakendId")
                .add("upperBreakendId")
                .add("lowerBreakendIsStart")
                .add("upperBreakendIsStart")
                .add("chromosome")
                .add("arm")
                .add("assembled")
                .add("traversedSVCount")
                .add("length")
                .add("ploidy")
                .add("pseudogeneInfo")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final LinxLink svData) {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(svData.clusterId()))
                .add(String.valueOf(svData.chainId()))
                .add(String.valueOf(svData.chainIndex()))
                .add(String.valueOf(svData.chainCount()))
                .add(String.valueOf(svData.lowerSvId()))
                .add(String.valueOf(svData.upperSvId()))
                .add(String.valueOf(svData.lowerBreakendIsStart()))
                .add(String.valueOf(svData.upperBreakendIsStart()))
                .add(String.valueOf(svData.chromosome()))
                .add(String.valueOf(svData.arm()))
                .add(String.valueOf(svData.assembled()))
                .add(String.valueOf(svData.traversedSVCount()))
                .add(String.valueOf(svData.length()))
                .add(String.valueOf(svData.ploidy()))
                .add(String.valueOf(svData.pseudogeneInfo()))
                .toString();
    }

    @NotNull
    private static LinxLink fromString(@NotNull final String tiData)
    {
        String[] values = tiData.split(DELIMITER);

        int index = 0;

        return ImmutableLinxLink.builder()
                .clusterId(Integer.valueOf(values[index++]))
                .chainId(Integer.valueOf(values[index++]))
                .chainIndex(Integer.valueOf(values[index++]))
                .chainCount(Integer.valueOf(values[index++]))
                .upperSvId(Integer.valueOf(values[index++]))
                .lowerSvId(Integer.valueOf(values[index++]))
                .lowerBreakendIsStart(Boolean.valueOf(values[index++]))
                .upperBreakendIsStart(Boolean.valueOf(values[index++]))
                .chromosome(values[index++])
                .arm(values[index++])
                .assembled(Boolean.valueOf(values[index++]))
                .traversedSVCount(Integer.valueOf(values[index++]))
                .length(Long.valueOf(values[index++]))
                .ploidy(Integer.valueOf(values[index++]))
                .pseudogeneInfo(values[index++])
                .build();
    }
}
