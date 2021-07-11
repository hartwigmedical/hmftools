package com.hartwig.hmftools.common.sv.linx;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.sv.linx.LinxCluster.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
public abstract class LinxBreakend
{
    public abstract int id();
    public abstract int svId();
    public abstract boolean isStart();
    public abstract String gene();
    public abstract String transcriptId();
    public abstract boolean canonical();
    public abstract String geneOrientation();
    public abstract boolean disruptive();
    public abstract boolean reportedDisruption();
    public abstract double undisruptedCopyNumber();
    public abstract String regionType();
    public abstract String codingContext();
    public abstract String biotype();
    public abstract int exonicBasePhase();
    public abstract int nextSpliceExonRank();
    public abstract int nextSpliceExonPhase();
    public abstract int nextSpliceDistance();
    public abstract int totalExonCount();

    // additional fields for patient report
    public abstract String type();
    public abstract String chromosome();
    public abstract int orientation();
    public abstract int strand();
    public abstract String chrBand();
    public abstract int exonUp();
    public abstract int exonDown();
    public abstract double junctionCopyNumber();

    private static final String FILE_EXTENSION = ".linx.breakend.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<LinxBreakend> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<LinxBreakend> breakends) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(breakends));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<LinxBreakend> breakends)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        breakends.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<LinxBreakend> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("id")).map(LinxBreakend::fromString).collect(toList());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("id")
                .add("svId")
                .add("isStart")
                .add("gene")
                .add("transcriptId")
                .add("canonical")
                .add("geneOrientation")
                .add("disruptive")
                .add("reportedDisruption")
                .add("undisruptedCopyNumber")
                .add("regionType")
                .add("codingContext")
                .add("biotype")
                .add("exonicBasePhase")
                .add("nextSpliceExonRank")
                .add("nextSpliceExonPhase")
                .add("nextSpliceDistance")
                .add("totalExonCount")
                .add("type")
                .add("chromosome")
                .add("orientation")
                .add("strand")
                .add("chrBand")
                .add("exonUp")
                .add("exonDown")
                .add("junctionCopyNumber")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final LinxBreakend breakend)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(breakend.id()))
                .add(String.valueOf(breakend.svId()))
                .add(String.valueOf(breakend.isStart()))
                .add(String.valueOf(breakend.gene()))
                .add(String.valueOf(breakend.transcriptId()))
                .add(String.valueOf(breakend.canonical()))
                .add(String.valueOf(breakend.geneOrientation()))
                .add(String.valueOf(breakend.disruptive()))
                .add(String.valueOf(breakend.reportedDisruption()))
                .add(String.format("%.4f", breakend.undisruptedCopyNumber()))
                .add(String.valueOf(breakend.regionType()))
                .add(String.valueOf(breakend.codingContext()))
                .add(String.valueOf(breakend.biotype()))
                .add(String.valueOf(breakend.exonicBasePhase()))
                .add(String.valueOf(breakend.nextSpliceExonRank()))
                .add(String.valueOf(breakend.nextSpliceExonPhase()))
                .add(String.valueOf(breakend.nextSpliceDistance()))
                .add(String.valueOf(breakend.totalExonCount()))
                .add(String.valueOf(breakend.type()))
                .add(String.valueOf(breakend.chromosome()))
                .add(String.valueOf(breakend.orientation()))
                .add(String.valueOf(breakend.strand()))
                .add(String.valueOf(breakend.chrBand()))
                .add(String.valueOf(breakend.exonUp()))
                .add(String.valueOf(breakend.exonDown()))
                .add(String.format("%.4f", breakend.junctionCopyNumber()))
                .toString();
    }

    @NotNull
    private static LinxBreakend fromString(@NotNull final String breakend)
    {
        String[] values = breakend.split(DELIMITER, -1);

        int index = 0;

        return ImmutableLinxBreakend.builder()
                .id(Integer.parseInt(values[index++]))
                .svId(Integer.parseInt(values[index++]))
                .isStart(Boolean.parseBoolean(values[index++]))
                .gene(values[index++])
                .transcriptId(values[index++])
                .canonical(Boolean.parseBoolean(values[index++]))
                .geneOrientation(values[index++])
                .disruptive(Boolean.parseBoolean(values[index++]))
                .reportedDisruption(Boolean.parseBoolean(values[index++]))
                .undisruptedCopyNumber(Double.parseDouble(values[index++]))
                .regionType(values[index++])
                .codingContext(values[index++])
                .biotype(values[index++])
                .exonicBasePhase(Integer.parseInt(values[index++]))
                .nextSpliceExonRank(Integer.parseInt(values[index++]))
                .nextSpliceExonPhase(Integer.parseInt(values[index++]))
                .nextSpliceDistance(Integer.parseInt(values[index++]))
                .totalExonCount(Integer.parseInt(values[index++]))
                .type(values[index++])
                .chromosome(values[index++])
                .orientation(Integer.parseInt(values[index++]))
                .strand(Integer.parseInt(values[index++]))
                .chrBand(values[index++])
                .exonUp(Integer.parseInt(values[index++]))
                .exonDown(Integer.parseInt(values[index++]))
                .junctionCopyNumber(Double.parseDouble(values[index++]))
                .build();
    }

}
