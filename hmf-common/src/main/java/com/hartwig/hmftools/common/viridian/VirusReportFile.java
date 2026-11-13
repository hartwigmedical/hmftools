package com.hartwig.hmftools.common.viridian;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.common.virus.VirusType;

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
                String reportingType = row.getStringOrNull(Column.interpretation);
                return new VirusReport(
                        row.get(Column.name),
                        row.get(Column.virus_group),
                        ReportQcStatus.valueOf(row.get(Column.qcStatus)),
                        reportingType != null ? VirusType.fromVirusName(reportingType) : null,
                        row.getInt(Column.integrations),
                        row.getDouble(Column.percentageCovered),
                        row.getDouble(Column.meanCoverage),
                        row.getDoubleOrNull(Column.copies_per_tumor_cell),
                        ClonalityStatus.valueOf(row.get(Column.clonality_status)),
                        row.getBoolean(Column.reported));
            }).toList();
        }
    }

    public static void write(String filename, List<VirusReport> reports)
    {
        DelimFileWriter.write(
                filename, Column.values(), reports, (report, row) ->
                {
                    row.set(Column.name, report.virusName());
                    row.set(Column.qcStatus, report.qcStatus().name());
                    row.set(Column.integrations, report.integrations());
                    row.setOrNull(Column.interpretation, report.reportingType() != null ? report.reportingType().toString() : null);
                    row.set(Column.percentageCovered, report.coveragePercentage());
                    row.set(Column.meanCoverage, report.meanDepth());
                    row.set(Column.reported, report.isReported());
                    row.set(Column.clonality_status, report.clonalityStatus().name());
                    row.setOrNull(Column.copies_per_tumor_cell, report.copiesPerTumorCell());
                    row.set(Column.virus_group, report.virusGroup());
                });
    }

    private enum Column
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
