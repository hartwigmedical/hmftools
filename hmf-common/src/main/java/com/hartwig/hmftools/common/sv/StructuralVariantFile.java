package com.hartwig.hmftools.common.sv;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    public static final String DELIMITER = "\t";
    private static final String FILE_EXTENSION = ".sv_data.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<StructuralVariantData> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<StructuralVariantData> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<StructuralVariantData> svDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        svDataList.stream().map(StructuralVariantFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<StructuralVariantData> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("id")).map(StructuralVariantFile::fromString).collect(toList());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("id")
                .add("startChromosome")
                .add("endChromosome")
                .add("startPosition")
                .add("endPosition")
                .add("startOrientation")
                .add("endOrientation")
                .add("startHomologySequence")
                .add("endHomologySequence")
                .add("startAF")
                .add("endAF")
                .add("junctionCopyNumber")
                .add("adjustedStartAF")
                .add("adjustedEndAF")
                .add("adjustedStartCopyNumber")
                .add("adjustedEndCopyNumber")
                .add("adjustedStartCopyNumberChange")
                .add("adjustedEndCopyNumberChange")
                .add("insertSequence")
                .add("type")
                .add("filter")
                .add("imprecise")
                .add("qualityScore")
                .add("event")
                .add("startTumorVariantFragmentCount")
                .add("startTumorReferenceFragmentCount")
                .add("startNormalVariantFragmentCount")
                .add("startNormalReferenceFragmentCount")
                .add("endTumorVariantFragmentCount")
                .add("endTumorReferenceFragmentCount")
                .add("endNormalVariantFragmentCount")
                .add("endNormalReferenceFragmentCount")
                .add("startIntervalOffsetStart")
                .add("startIntervalOffsetEnd")
                .add("endIntervalOffsetStart")
                .add("endIntervalOffsetEnd")
                .add("inexactHomologyOffsetStart")
                .add("inexactHomologyOffsetEnd")
                .add("startLinkedBy")
                .add("endLinkedBy")
                .add("vcfId")
                .add("recovered")
                .add("recoveryMethod")
                .add("recoveryFilter")
                .add("startRefContext")
                .add("endRefContext")
                .add("insertSequenceAlignments")
                .add("insertSequenceRepeatClass")
                .add("insertSequenceRepeatType")
                .add("insertSequenceRepeatOrientation")
                .add("insertSequenceRepeatCoverage")
                .add("startAnchoringSupportDistance")
                .add("endAnchoringSupportDistance")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final StructuralVariantData svData)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(svData.id()))
                .add(String.valueOf(svData.startChromosome()))
                .add(String.valueOf(svData.endChromosome()))
                .add(String.valueOf(svData.startPosition()))
                .add(String.valueOf(svData.endPosition()))
                .add(String.valueOf(svData.startOrientation()))
                .add(String.valueOf(svData.endOrientation()))
                .add(String.valueOf(svData.startHomologySequence()))
                .add(String.valueOf(svData.endHomologySequence()))
                .add(FORMAT.format(svData.startAF()))
                .add(FORMAT.format(svData.endAF()))
                .add(FORMAT.format(svData.junctionCopyNumber()))
                .add(FORMAT.format(svData.adjustedStartAF()))
                .add(FORMAT.format(svData.adjustedEndAF()))
                .add(FORMAT.format(svData.adjustedStartCopyNumber()))
                .add(FORMAT.format(svData.adjustedEndCopyNumber()))
                .add(FORMAT.format(svData.adjustedStartCopyNumberChange()))
                .add(FORMAT.format(svData.adjustedEndCopyNumberChange()))
                .add(String.valueOf(svData.insertSequence()))
                .add(String.valueOf(svData.type()))
                .add(String.valueOf(svData.filter()))
                .add(String.valueOf(svData.imprecise()))
                .add(FORMAT.format(svData.qualityScore()))
                .add(String.valueOf(svData.event()))
                .add(String.valueOf(svData.startTumorVariantFragmentCount()))
                .add(String.valueOf(svData.startTumorReferenceFragmentCount()))
                .add(String.valueOf(svData.startNormalVariantFragmentCount()))
                .add(String.valueOf(svData.startNormalReferenceFragmentCount()))
                .add(String.valueOf(svData.endTumorVariantFragmentCount()))
                .add(String.valueOf(svData.endTumorReferenceFragmentCount()))
                .add(String.valueOf(svData.endNormalVariantFragmentCount()))
                .add(String.valueOf(svData.endNormalReferenceFragmentCount()))
                .add(String.valueOf(svData.startIntervalOffsetStart()))
                .add(String.valueOf(svData.startIntervalOffsetEnd()))
                .add(String.valueOf(svData.endIntervalOffsetStart()))
                .add(String.valueOf(svData.endIntervalOffsetEnd()))
                .add(String.valueOf(svData.inexactHomologyOffsetStart()))
                .add(String.valueOf(svData.inexactHomologyOffsetEnd()))
                .add(String.valueOf(svData.startLinkedBy()))
                .add(String.valueOf(svData.endLinkedBy()))
                .add(String.valueOf(svData.vcfId()))
                .add(String.valueOf(svData.recovered()))
                .add(String.valueOf(svData.recoveryMethod()))
                .add(String.valueOf(svData.recoveryFilter()))
                .add(String.valueOf(svData.startRefContext()))
                .add(String.valueOf(svData.endRefContext()))
                .add(String.valueOf(svData.insertSequenceAlignments()))
                .add(String.valueOf(svData.insertSequenceRepeatClass()))
                .add(String.valueOf(svData.insertSequenceRepeatType()))
                .add(String.valueOf(svData.insertSequenceRepeatOrientation()))
                .add(String.valueOf(svData.insertSequenceRepeatCoverage()))
                .add(String.valueOf(svData.startAnchoringSupportDistance()))
                .add(String.valueOf(svData.endAnchoringSupportDistance()))
                .toString();
    }

    @NotNull
    public static StructuralVariantData fromString(@NotNull final String svData)
    {
        String[] values = svData.split(DELIMITER);

        int index = 0;

        final ImmutableStructuralVariantData.Builder builder = ImmutableStructuralVariantData.builder()
                .id(Integer.parseInt(values[index++]))
                .startChromosome(values[index++])
                .endChromosome(getStrValue(values[index++]))
                .startPosition(Integer.parseInt(values[index++]))
                .endPosition(getIntValue(values[index++]))
                .startOrientation(Byte.parseByte(values[index++]))
                .endOrientation(getByteValue(values[index++]))
                .startHomologySequence(values[index++])
                .endHomologySequence(getStrValue(values[index++]))
                .startAF(Double.parseDouble(values[index++]))
                .endAF(getDoubleValue(values[index++]))
                .junctionCopyNumber(Double.parseDouble(values[index++]))
                .adjustedStartAF(Double.parseDouble(values[index++]))
                .adjustedEndAF(getDoubleValue(values[index++]))
                .adjustedStartCopyNumber(Double.parseDouble(values[index++]))
                .adjustedEndCopyNumber(getDoubleValue(values[index++]))
                .adjustedStartCopyNumberChange(Double.parseDouble(values[index++]))
                .adjustedEndCopyNumberChange(getDoubleValue(values[index++]))
                .insertSequence(values[index++]);

        StructuralVariantType type = StructuralVariantType.fromAttribute(values[index++]);
        final String filterStr = values[index++];
        if(type == SGL && filterStr.equals(INFERRED))
            type = INF;

        builder.type(type)
                .filter(filterStr)
                .imprecise(Boolean.parseBoolean(values[index++]))
                .qualityScore(Double.parseDouble(values[index++]))
                .event(values[index++])
                .startTumorVariantFragmentCount(Integer.parseInt(values[index++]))
                .startTumorReferenceFragmentCount(Integer.parseInt(values[index++]))
                .startNormalVariantFragmentCount(Integer.parseInt(values[index++]))
                .startNormalReferenceFragmentCount(Integer.parseInt(values[index++]))
                .endTumorVariantFragmentCount(getIntValue(values[index++]))
                .endTumorReferenceFragmentCount(getIntValue(values[index++]))
                .endNormalVariantFragmentCount(getIntValue(values[index++]))
                .endNormalReferenceFragmentCount(getIntValue(values[index++]))
                .startIntervalOffsetStart(Integer.parseInt(values[index++]))
                .startIntervalOffsetEnd(Integer.parseInt(values[index++]))
                .endIntervalOffsetStart(getIntValue(values[index++]))
                .endIntervalOffsetEnd(getIntValue(values[index++]))
                .inexactHomologyOffsetStart(Integer.parseInt(values[index++]))
                .inexactHomologyOffsetEnd(getIntValue(values[index++]))
                .startLinkedBy(values[index++])
                .endLinkedBy(values[index++])
                .vcfId(values[index++])
                .startRefContext(values[index++])
                .endRefContext(values[index++])
                .recovered(Boolean.parseBoolean(values[index++]))
                .recoveryMethod(values[index++])
                .recoveryFilter(values[index++])
                .insertSequenceAlignments(values[index++])
                .insertSequenceRepeatClass(values[index++])
                .insertSequenceRepeatType(values[index++])
                .insertSequenceRepeatOrientation(getByteValue(values[index++]))
                .insertSequenceRepeatCoverage(getDoubleValue(values[index++]))
                .startAnchoringSupportDistance(Integer.parseInt(values[index++]))
                .endAnchoringSupportDistance(getIntValue(values[index++]));

        return builder.build();
    }

    private static String getStrValue(final String value) { return value.equals("NULL") ? "" : value; }
    private static double getDoubleValue(final String value) { return value.equals("NULL") ? 0 : Double.parseDouble(value); }
    private static int getIntValue(final String value) { return value.equals("NULL") ? 0 : Integer.parseInt(value); }
    private static long getLongValue(final String value) { return value.equals("NULL") ? 0 : Long.parseLong(value); }
    private static byte getByteValue(final String value) { return value.equals("NULL") ? 0 : Byte.valueOf(value); }
}
