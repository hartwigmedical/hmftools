package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.lims.LimsModel;
import com.hartwig.hmftools.patientreporter.report.ReportWriter;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.exception.DRException;

public class NotSequenceableReporter {

    @NotNull
    private final CpctEcrfModel cpctEcrfModel;
    @NotNull
    private final LimsModel limsModel;
    @NotNull
    private final ReportWriter reportWriter;

    public NotSequenceableReporter(@NotNull final CpctEcrfModel cpctEcrfModel, @NotNull final LimsModel limsModel,
            @NotNull final ReportWriter reportWriter) {
        this.cpctEcrfModel = cpctEcrfModel;
        this.limsModel = limsModel;
        this.reportWriter = reportWriter;
    }

    public void run(@NotNull final String sample, @NotNull final NotSequenceableReason reason, @NotNull final NotSequenceableStudy study)
            throws IOException, DRException {
        final String tumorType = PatientReporterHelper.extractTumorType(cpctEcrfModel, sample);
        final String tumorPercentageString = PatientReportFormat.formatNullablePercent(limsModel.findTumorPercentageForSample(sample));
        reportWriter.writeNonSequenceableReport(sample, tumorType, tumorPercentageString, reason, study);
    }
}
