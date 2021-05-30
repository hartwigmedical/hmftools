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
public abstract class LinxLink
{
    public abstract int clusterId();
    public abstract int chainId();
    public abstract String chainIndex();
    public abstract int chainCount();
    public abstract int lowerSvId();
    public abstract int upperSvId();
    public abstract boolean lowerBreakendIsStart();
    public abstract boolean upperBreakendIsStart();
    public abstract String chromosome();
    public abstract String arm();
    public abstract boolean assembled();
    public abstract int traversedSVCount();
    public abstract long length();
    public abstract double junctionCopyNumber();
    public abstract double junctionCopyNumberUncertainty();
    public abstract String pseudogeneInfo();
    public abstract boolean ecDna();

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
    private static List<String> toLines(@NotNull final List<LinxLink> svDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        svDataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<LinxLink> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("clusterId")).map(LinxLink::fromString).collect(toList());
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
                .add("junctionCopyNumber")
                .add("junctionCopyNumberUncertainty")
                .add("pseudogeneInfo")
                .add("ecDna")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final LinxLink svData)
    {
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
                .add(String.valueOf(svData.junctionCopyNumber()))
                .add(String.valueOf(svData.junctionCopyNumberUncertainty()))
                .add(String.valueOf(svData.pseudogeneInfo()))
                .add(String.valueOf(svData.ecDna()))
                .toString();
    }

    @NotNull
    private static LinxLink fromString(@NotNull final String tiData)
    {
        String[] values = tiData.split(DELIMITER);

        int index = 0;

        return ImmutableLinxLink.builder()
                .clusterId(Integer.parseInt(values[index++]))
                .chainId(Integer.parseInt(values[index++]))
                .chainIndex(values[index++])
                .chainCount(Integer.parseInt(values[index++]))
                .lowerSvId(Integer.parseInt(values[index++]))
                .upperSvId(Integer.parseInt(values[index++]))
                .lowerBreakendIsStart(Boolean.parseBoolean(values[index++]))
                .upperBreakendIsStart(Boolean.parseBoolean(values[index++]))
                .chromosome(values[index++])
                .arm(values[index++])
                .assembled(Boolean.parseBoolean(values[index++]))
                .traversedSVCount(Integer.parseInt(values[index++]))
                .length(Long.parseLong(values[index++]))
                .junctionCopyNumber(Double.parseDouble(values[index++]))
                .junctionCopyNumberUncertainty(Double.parseDouble(values[index++]))
                .pseudogeneInfo(index < values.length ? values[index++] : "")
                .ecDna(Boolean.parseBoolean(values[index++]))
                .build();
    }

}
