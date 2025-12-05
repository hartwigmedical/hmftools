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
import com.hartwig.hmftools.common.driver.DriverType;

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
                .add("relativeMinCopyNumber")
                .add("depthWindowCount")
                .add("gcContent")
                .add("reportedStatus")
                .add("driverType")
                .toString();
    }

    private static String toString(final GeneCopyNumber geneCopyNumber)
    {
        return new StringJoiner(TSV_DELIM).add(geneCopyNumber.chromosome())
                .add(String.valueOf(geneCopyNumber.start()))
                .add(String.valueOf(geneCopyNumber.end()))
                .add(geneCopyNumber.geneName())
                .add(FORMAT.format(geneCopyNumber.MinCopyNumber))
                .add(FORMAT.format(geneCopyNumber.MaxCopyNumber))
                .add(String.valueOf(geneCopyNumber.SomaticRegions))
                .add(geneCopyNumber.TransName)
                .add(String.valueOf(geneCopyNumber.IsCanonical))
                .add(geneCopyNumber.ChromosomeBand)
                .add(String.valueOf(geneCopyNumber.MinRegions))
                .add(String.valueOf(geneCopyNumber.MinRegionStart))
                .add(String.valueOf(geneCopyNumber.MinRegionEnd))
                .add(String.valueOf(geneCopyNumber.MinRegionStartSupport))
                .add(String.valueOf(geneCopyNumber.MinRegionEndSupport))
                .add(String.valueOf(geneCopyNumber.MinRegionMethod))
                .add(FORMAT.format(geneCopyNumber.MinMinorAlleleCopyNumber))
                .add(FORMAT.format(geneCopyNumber.RelativeMinCopyNumber))
                .add(String.valueOf(geneCopyNumber.DepthWindowCount))
                .add(FORMAT.format(geneCopyNumber.GcContent))
                .add(String.valueOf(geneCopyNumber.reportedStatus()))
                .add(String.valueOf(geneCopyNumber.driverType()))
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
        int canonicalIndex = fieldsIndexMap.get("isCanonical");
        int chrBandIndex = fieldsIndexMap.get("chromosomeBand");
        int minRegionIndex = fieldsIndexMap.get("minRegions");
        int minRegionStartIndex = fieldsIndexMap.get("minRegionStart");
        int minRegionEndIndex = fieldsIndexMap.get("minRegionEnd");
        int minRegionStartSupIndex = fieldsIndexMap.get("minRegionStartSupport");
        int minRegionEndSupIndex = fieldsIndexMap.get("minRegionEndSupport");
        int minRegionMethodIndex = fieldsIndexMap.get("minRegionMethod");
        int mmACnIndex = fieldsIndexMap.get("minMinorAlleleCopyNumber");
        int dwcIndex = fieldsIndexMap.get("depthWindowCount");
        Integer gcIndex = fieldsIndexMap.get("gcContent");
        Integer reportedIndex = fieldsIndexMap.get("reportedStatus");
        Integer driverIndex = fieldsIndexMap.get("driverType");
        Integer relMinCnIndex = fieldsIndexMap.get("relativeMinCopyNumber");

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();

        for(final String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            double gcContent = gcIndex != null ? Double.parseDouble(values[gcIndex]) : 0; // addded in Purple v4.2
            double relativeMinCopyNumber = relMinCnIndex != null ? Double.parseDouble(values[relMinCnIndex]) : 0; // added in Purple v4.3

            GeneCopyNumber geneCopyNumber = new GeneCopyNumber(
                    values[chrIndex], Integer.parseInt(values[startIndex]), Integer.parseInt(values[endIndex]),
                    values[geneIndex], values[transIdIndex], Boolean.parseBoolean(values[canonicalIndex]), values[chrBandIndex],
                    Double.parseDouble(values[maxCnIndex]), Double.parseDouble(values[minCnIndex]), Double.parseDouble(values[mmACnIndex]),
                    Integer.parseInt(values[somRegionsIndex]), Integer.parseInt(values[minRegionIndex]),
                    Integer.parseInt(values[minRegionStartIndex]), Integer.parseInt(values[minRegionEndIndex]),
                    Integer.parseInt(values[dwcIndex]), gcContent,
                    SegmentSupport.valueOf(values[minRegionStartSupIndex]), SegmentSupport.valueOf(values[minRegionEndSupIndex]),
                    CopyNumberMethod.valueOf(values[minRegionMethodIndex]), relativeMinCopyNumber);

            // added in Purple v4.3
            if(reportedIndex != null)
            {
                geneCopyNumber.setReportedStatus(ReportedStatus.valueOf(values[reportedIndex]));
            }

            if(driverIndex != null)
            {
                geneCopyNumber.setDriverType(DriverType.valueOf(values[driverIndex]));
            }

            geneCopyNumbers.add(geneCopyNumber);
        }

        return geneCopyNumbers;
    }
}
