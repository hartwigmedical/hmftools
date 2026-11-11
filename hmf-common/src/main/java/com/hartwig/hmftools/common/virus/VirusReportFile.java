package com.hartwig.hmftools.common.virus;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.jetbrains.annotations.Nullable;

public final class VirusReportFile
{
    private static final String FILE_EXTENSION = ".viridian.report.tsv";

    public static String generateFilename(String basePath, String sampleId)
    {
        return generateFilename(basePath, sampleId, null);
    }

    public static String generateFilename(String basePath, String sampleId, @Nullable String outputId)
    {
        String filename = checkAddDirSeparator(basePath) + sampleId;
        if(outputId != null)
        {
            filename += "." + outputId;
        }
        return filename + FILE_EXTENSION;
    }

    public static List<VirusReport> read(String filename)
    {
        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            return reader.stream().map(row ->
            {
                String reportingType = row.getStringOrNull(Columns.interpretation);
                return new VirusReport(
                        row.get(Columns.name),
                        row.get(Columns.virus_group),
                        ReportQcStatus.valueOf(row.get(Columns.qcStatus)),
                        reportingType != null ? VirusType.fromVirusName(reportingType) : null,
                        row.getInt(Columns.integrations),
                        row.getDouble(Columns.percentageCovered),
                        row.getDouble(Columns.meanCoverage),
                        row.getDoubleOrNull(Columns.copies_per_tumor_cell),
                        ClonalityStatus.valueOf(row.get(Columns.clonality_status)),
                        row.getBoolean(Columns.reported));
            }).toList();
        }
    }

    public static void write(String filename, List<VirusReport> reports)
    {
        DelimFileWriter.write(
                filename, Columns.values(), reports, (report, row) ->
                {
                    row.set(Columns.name, report.virusName());
                    row.set(Columns.qcStatus, report.qcStatus().name());
                    row.set(Columns.integrations, report.integrations());
                    row.setOrNull(Columns.interpretation, report.reportingType() != null ? report.reportingType().toString() : null);
                    row.set(Columns.percentageCovered, report.coveragePercentage());
                    row.set(Columns.meanCoverage, report.meanDepth());
                    row.set(Columns.reported, report.isReported());
                    row.set(Columns.clonality_status, report.clonalityStatus().name());
                    row.setOrNull(Columns.copies_per_tumor_cell, report.copiesPerTumorCell());
                    row.set(Columns.virus_group, report.virusGroup());
                });
    }

    private enum Columns
    {
        name,
        qcStatus,
        integrations,
        interpretation,
        percentageCovered,
        meanCoverage,
        reported,
        clonality_status,
        copies_per_tumor_cell,
        virus_group,
    }
}
