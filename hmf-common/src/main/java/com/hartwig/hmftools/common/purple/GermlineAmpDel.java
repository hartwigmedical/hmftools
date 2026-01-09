package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class GermlineAmpDel
{
    public final String GeneName;
    public final String Chromosome;
    public final String ChromosomeBand;
    public final int RegionStart;
    public final int RegionEnd;
    public final int DepthWindowCount;
    public final int ExonStart;
    public final int ExonEnd;
    public final GermlineDetectionMethod DetectionMethod;
    public final GermlineStatus NormalStatus;
    public final GermlineStatus TumorStatus;
    public final double GermlineCopyNumber;
    public final double TumorCopyNumber;
    public final String Filter;
    public final int CohortFrequency;
    public final ReportedStatus Reported;

    public GermlineAmpDel(
            final String geneName, final String chromosome, final String chromosomeBand, final int regionStart, final int regionEnd,
            final int depthWindowCount, final int exonStart, final int exonEnd, final GermlineDetectionMethod detectionMethod,
            final GermlineStatus normalStatus, final GermlineStatus tumorStatus, final double germlineCopyNumber, final double tumorCopyNumber,
            final String filter, final int cohortFrequency, final ReportedStatus reportedStatus)
    {
        GeneName = geneName;
        Chromosome = chromosome;
        ChromosomeBand = chromosomeBand;
        RegionStart = regionStart;
        RegionEnd = regionEnd;
        DepthWindowCount = depthWindowCount;
        ExonStart = exonStart;
        ExonEnd = exonEnd;
        DetectionMethod = detectionMethod;
        NormalStatus = normalStatus;
        TumorStatus = tumorStatus;
        TumorCopyNumber = tumorCopyNumber;
        GermlineCopyNumber = germlineCopyNumber;
        Filter = filter;
        CohortFrequency = cohortFrequency;
        Reported = reportedStatus;
    }

    private static final String EXTENSION = ".purple.germline.deletion.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static List<GermlineAmpDel> read(final String fileName) throws IOException
    {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(final String fileName, List<GermlineAmpDel> deletions) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(deletions));
    }

    private static List<String> toLines(final List<GermlineAmpDel> deletions)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        deletions.forEach(x -> lines.add(toString(x)));
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM, "", "")
                .add("gene")
                .add("chromosome")
                .add("chromosomeBand")
                .add("regionStart")
                .add("regionEnd")
                .add("depthWindowCount")
                .add("exonStart")
                .add("exonEnd")
                .add("detectionMethod")
                .add("germlineStatus")
                .add("tumorStatus")
                .add("germlineCopyNumber")
                .add("tumorCopyNumber")
                .add("filter")
                .add("cohortFrequency")
                .add("reportedStatus")
                .toString();
    }

    private static String toString(final GermlineAmpDel deletion)
    {
        return new StringJoiner(TSV_DELIM)
                .add(deletion.GeneName)
                .add(deletion.Chromosome)
                .add(deletion.ChromosomeBand)
                .add(String.valueOf(deletion.RegionStart))
                .add(String.valueOf(deletion.RegionEnd))
                .add(String.valueOf(deletion.DepthWindowCount))
                .add(String.valueOf(deletion.ExonStart))
                .add(String.valueOf(deletion.ExonEnd))
                .add(deletion.DetectionMethod.toString())
                .add(deletion.NormalStatus.toString())
                .add(deletion.TumorStatus.toString())
                .add(String.format("%.2f", deletion.GermlineCopyNumber))
                .add(String.format("%.2f", deletion.TumorCopyNumber))
                .add(deletion.Filter)
                .add(String.valueOf(deletion.CohortFrequency))
                .add(String.valueOf(deletion.Reported))
                .toString();
    }

    static List<GermlineAmpDel> fromLines(final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        List<GermlineAmpDel> deletions = Lists.newArrayList();

        Integer reportedStatusIndex = fieldsIndexMap.get("reportedStatus");
        Integer reportedIndex = fieldsIndexMap.get("reported");

        for(final String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            ReportedStatus reportedStatus;

            if(reportedStatusIndex != null)
            {
                reportedStatus = ReportedStatus.valueOf(values[reportedStatusIndex]);
            }
            else
            {
                reportedStatus = Boolean.parseBoolean(values[reportedIndex]) ? ReportedStatus.REPORTED : ReportedStatus.NONE;
            }

            deletions.add(new GermlineAmpDel(
                    values[fieldsIndexMap.get("gene")],
                    values[fieldsIndexMap.get("chromosome")],
                    fieldsIndexMap.containsKey("chromosomeBand") ? values[fieldsIndexMap.get("chromosomeBand")] : "",
                    Integer.parseInt(values[fieldsIndexMap.get("regionStart")]), Integer.parseInt(values[fieldsIndexMap.get("regionEnd")]),
                    Integer.parseInt(values[fieldsIndexMap.get("depthWindowCount")]),
                    Integer.parseInt(values[fieldsIndexMap.get("exonStart")]), Integer.parseInt(values[fieldsIndexMap.get("exonEnd")]),
                    GermlineDetectionMethod.valueOf(values[fieldsIndexMap.get("detectionMethod")]),
                    GermlineStatus.valueOf(values[fieldsIndexMap.get("germlineStatus")]),
                    GermlineStatus.valueOf(values[fieldsIndexMap.get("tumorStatus")]),
                    Double.parseDouble(values[fieldsIndexMap.get("germlineCopyNumber")]),
                    Double.parseDouble(values[fieldsIndexMap.get("tumorCopyNumber")]),
                    values[fieldsIndexMap.get("filter")],
                    Integer.parseInt(values[fieldsIndexMap.get("cohortFrequency")]),
                    reportedStatus));
        }

        return deletions;
    }
}
