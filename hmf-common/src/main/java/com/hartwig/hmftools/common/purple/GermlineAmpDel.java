package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class GermlineAmpDel
{
    public final String GeneName;
    public final String Transcript;
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
            final String geneName, final String transcript, final String chromosome, final String chromosomeBand, int regionStart, int regionEnd,
            final int depthWindowCount, final int exonStart, final int exonEnd, final GermlineDetectionMethod detectionMethod,
            final GermlineStatus normalStatus, final GermlineStatus tumorStatus, final double germlineCopyNumber, final double tumorCopyNumber,
            final String filter, final int cohortFrequency, final ReportedStatus reportedStatus)
    {
        GeneName = geneName;
        Transcript = transcript;
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

    private static final String EXTENSION = ".purple.germline_amp_del.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static List<GermlineAmpDel> read(final String fileName) throws IOException
    {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(final String fileName, List<GermlineAmpDel> ampDels) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(ampDels));
    }

    private static List<String> toLines(final List<GermlineAmpDel> ampDels)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        ampDels.forEach(x -> lines.add(toString(x)));
        return lines;
    }

    private enum Columns
    {
        gene,
        transcript,
        chromosome,
        chromosomeBand,
        regionStart,
        regionEnd,
        depthWindowCount,
        exonStart,
        exonEnd,
        detectionMethod,
        germlineStatus,
        tumorStatus,
        germlineCopyNumber,
        tumorCopyNumber,
        filter,
        cohortFrequency,
        reportedStatus;
    }

    private static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        Arrays.stream(Columns.values()).forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    private static String toString(final GermlineAmpDel ampDel)
    {
        return new StringJoiner(TSV_DELIM)
                .add(ampDel.GeneName)
                .add(ampDel.Transcript)
                .add(ampDel.Chromosome)
                .add(ampDel.ChromosomeBand)
                .add(String.valueOf(ampDel.RegionStart))
                .add(String.valueOf(ampDel.RegionEnd))
                .add(String.valueOf(ampDel.DepthWindowCount))
                .add(String.valueOf(ampDel.ExonStart))
                .add(String.valueOf(ampDel.ExonEnd))
                .add(ampDel.DetectionMethod.toString())
                .add(ampDel.NormalStatus.toString())
                .add(ampDel.TumorStatus.toString())
                .add(String.format("%.2f", ampDel.GermlineCopyNumber))
                .add(String.format("%.2f", ampDel.TumorCopyNumber))
                .add(ampDel.Filter)
                .add(String.valueOf(ampDel.CohortFrequency))
                .add(String.valueOf(ampDel.Reported))
                .toString();
    }

    static List<GermlineAmpDel> fromLines(final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        lines.remove(0);

        List<GermlineAmpDel> ampDels = Lists.newArrayList();

        Integer reportedStatusIndex = fieldsIndexMap.get(Columns.reportedStatus.toString());

        // for pre-v3.0
        Integer reportedIndex = fieldsIndexMap.get("reported");
        Integer transcriptIndex = fieldsIndexMap.get(Columns.transcript.toString());

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

            ampDels.add(new GermlineAmpDel(
                    values[fieldsIndexMap.get(Columns.gene.toString())],
                    transcriptIndex != null ? values[fieldsIndexMap.get(Columns.transcript.toString())] : "",
                    values[fieldsIndexMap.get(Columns.chromosome.toString())],
                    values[fieldsIndexMap.get(Columns.chromosomeBand.toString())],
                    Integer.parseInt(values[fieldsIndexMap.get(Columns.regionStart.toString())]),
                    Integer.parseInt(values[fieldsIndexMap.get(Columns.regionEnd.toString())]),
                    Integer.parseInt(values[fieldsIndexMap.get(Columns.depthWindowCount.toString())]),
                    Integer.parseInt(values[fieldsIndexMap.get(Columns.exonStart.toString())]),
                    Integer.parseInt(values[fieldsIndexMap.get(Columns.exonEnd.toString())]),
                    GermlineDetectionMethod.valueOf(values[fieldsIndexMap.get(Columns.detectionMethod.toString())]),
                    GermlineStatus.valueOf(values[fieldsIndexMap.get(Columns.germlineStatus.toString())]),
                    GermlineStatus.valueOf(values[fieldsIndexMap.get(Columns.tumorStatus.toString())]),
                    Double.parseDouble(values[fieldsIndexMap.get(Columns.germlineCopyNumber.toString())]),
                    Double.parseDouble(values[fieldsIndexMap.get(Columns.tumorCopyNumber.toString())]),
                    values[fieldsIndexMap.get(Columns.filter.toString())],
                    Integer.parseInt(values[fieldsIndexMap.get(Columns.cohortFrequency.toString())]),
                    reportedStatus));
        }

        return ampDels;
    }
}
