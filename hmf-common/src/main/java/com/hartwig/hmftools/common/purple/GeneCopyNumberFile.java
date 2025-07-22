package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

public final class GeneCopyNumberFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private static final String EXTENSION = ".purple.cnv.gene.tsv";

    private GeneCopyNumberFile() {}

    public static String generateFilenameForWriting(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    public static List<GeneCopyNumber> read(final String fileName) throws IOException
    {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(final String fileName, List<GeneCopyNumber> geneCopyNumbers) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(geneCopyNumbers));
    }

    @VisibleForTesting
    static List<String> toLines(final List<GeneCopyNumber> ratio)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        ratio.stream().map(GeneCopyNumberFile::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM, "", "")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("gene")
                .add("minCopyNumber")
                .add("maxCopyNumber")
                .add("somaticRegions")
                .add("transcriptId")
                .add("isCanonical")
                .add("chromosomeBand")
                .add("minRegions")
                .add("minRegionStart")
                .add("minRegionEnd")
                .add("minRegionStartSupport")
                .add("minRegionEndSupport")
                .add("minRegionMethod")
                .add("minMinorAlleleCopyNumber")
                .add("depthWindowCount")
                .add("gcContent")
                .toString();
    }

    private static String toString(final GeneCopyNumber geneCopyNumber)
    {
        return new StringJoiner(TSV_DELIM).add(geneCopyNumber.chromosome())
                .add(String.valueOf(geneCopyNumber.start()))
                .add(String.valueOf(geneCopyNumber.end()))
                .add(geneCopyNumber.geneName())
                .add(FORMAT.format(geneCopyNumber.minCopyNumber()))
                .add(FORMAT.format(geneCopyNumber.maxCopyNumber()))
                .add(String.valueOf(geneCopyNumber.somaticRegions()))
                .add(geneCopyNumber.transName())
                .add(String.valueOf(geneCopyNumber.isCanonical()))
                .add(geneCopyNumber.chromosomeBand())
                .add(String.valueOf(geneCopyNumber.minRegions()))
                .add(String.valueOf(geneCopyNumber.minRegionStart()))
                .add(String.valueOf(geneCopyNumber.minRegionEnd()))
                .add(String.valueOf(geneCopyNumber.minRegionStartSupport()))
                .add(String.valueOf(geneCopyNumber.minRegionEndSupport()))
                .add(String.valueOf(geneCopyNumber.minRegionMethod()))
                .add(FORMAT.format(geneCopyNumber.minMinorAlleleCopyNumber()))
                .add(String.valueOf(geneCopyNumber.depthWindowCount()))
                .add(FORMAT.format(geneCopyNumber.gcContent()))
                .toString();
    }

    @VisibleForTesting
    static List<GeneCopyNumber> fromLines(final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        int chrIndex = fieldsIndexMap.get("chromosome");
        int startIndex = fieldsIndexMap.get("start");
        int endIndex = fieldsIndexMap.get("end");
        int geneIndex = fieldsIndexMap.get("gene");
        int minCnIndex = fieldsIndexMap.get("minCopyNumber");
        int maxCnIndex = fieldsIndexMap.get("maxCopyNumber");
        int somRegionsIndex = fieldsIndexMap.get("somaticRegions");
        int transIdIndex = fieldsIndexMap.get("transcriptId");
        Integer canonicalIndex = fieldsIndexMap.get("isCanonical");
        int chrBandIndex = fieldsIndexMap.get("chromosomeBand");
        int minRegionIndex = fieldsIndexMap.get("minRegions");
        int minRegionStartIndex = fieldsIndexMap.get("minRegionStart");
        int minRegionEndIndex = fieldsIndexMap.get("minRegionEnd");
        int minRegionStartSupIndex = fieldsIndexMap.get("minRegionStartSupport");
        int minRegionEndSupIndex = fieldsIndexMap.get("minRegionEndSupport");
        int minRegionMethodIndex = fieldsIndexMap.get("minRegionMethod");
        int mmACnIndex = fieldsIndexMap.get("minMinorAlleleCopyNumber");
        Integer dwcIndex = fieldsIndexMap.get("depthWindowCount");
        Integer gcIndex = fieldsIndexMap.get("gcContent");

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();

        for(final String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            final ImmutableGeneCopyNumber.Builder builder = ImmutableGeneCopyNumber.builder()
                    .chromosome(values[chrIndex])
                    .start(Integer.parseInt(values[startIndex]))
                    .end(Integer.parseInt(values[endIndex]))
                    .geneName(values[geneIndex])
                    .minCopyNumber(Double.parseDouble(values[minCnIndex]))
                    .maxCopyNumber(Double.parseDouble(values[maxCnIndex]))
                    .somaticRegions(Integer.parseInt(values[somRegionsIndex]))
                    .transName(values[transIdIndex])
                    .isCanonical(canonicalIndex != null ? Boolean.parseBoolean(values[canonicalIndex]) : true)
                    .chromosomeBand(values[chrBandIndex])
                    .minRegions(minRegionIndex)
                    .minRegionStart(Integer.parseInt(values[minRegionStartIndex]))
                    .minRegionEnd(Integer.parseInt(values[minRegionEndIndex]))
                    .minRegionMethod(CopyNumberMethod.valueOf(values[minRegionMethodIndex]))
                    .minRegionStartSupport(SegmentSupport.valueOf(values[minRegionStartSupIndex]))
                    .minRegionEndSupport(SegmentSupport.valueOf(values[minRegionEndSupIndex]))
                    .minMinorAlleleCopyNumber(Double.parseDouble(values[mmACnIndex]))
                    .depthWindowCount(dwcIndex != null ? Integer.parseInt(values[dwcIndex]) : 0)
                    .gcContent(gcIndex != null ? Double.parseDouble(values[gcIndex]) : -1.0);

            geneCopyNumbers.add(builder.build());
        }

        return geneCopyNumbers;
    }
}
