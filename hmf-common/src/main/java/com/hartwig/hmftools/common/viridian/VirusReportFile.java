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
                        row.get(Column.virusGroup),
                        ReportQcStatus.valueOf(row.get(Column.qcStatus)),
                        reportingType != null ? VirusType.fromVirusName(reportingType) : null,
                        row.getInt(Column.integrations),
                        row.getDouble(Column.percentageCovered),
                        row.getDouble(Column.meanCoverage),
                        row.getDoubleOrNull(Column.copiesPerTumorCell),
                        ClonalityStatus.valueOf(row.get(Column.clonalityStatus)),
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
                    row.set(Column.clonalityStatus, report.clonalityStatus().name());
                    row.setOrNull(Column.copiesPerTumorCell, report.copiesPerTumorCell());
                    row.set(Column.virusGroup, report.virusGroup());
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
        clonalityStatus,
        copiesPerTumorCell,
        virusGroup,
    }
}
