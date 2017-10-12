package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.util.Optional;

import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.lims.LimsJsonModel;
import com.hartwig.hmftools.patientreporter.ImmutableNotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.NotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.report.ReportWriter;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.exception.DRException;

public class NotSequenceableReporter {

    @NotNull
    private final CpctEcrfModel cpctEcrfModel;
    @NotNull
    private final LimsJsonModel limsModel;
    @NotNull
    private final ReportWriter reportWriter;
    @NotNull
    private final CenterModel centerModel;

    public NotSequenceableReporter(@NotNull final CpctEcrfModel cpctEcrfModel, @NotNull final LimsJsonModel limsModel,
            @NotNull final CenterModel centerModel, @NotNull final ReportWriter reportWriter) {
        this.cpctEcrfModel = cpctEcrfModel;
        this.limsModel = limsModel;
        this.centerModel = centerModel;
        this.reportWriter = reportWriter;
    }

    public void run(@NotNull final String sample, @NotNull final NotSequenceableReason reason, @NotNull final NotSequenceableStudy study,
            @NotNull final String signaturePath, @Nullable final String comments) throws IOException, DRException {
        final String tumorType = PatientReporterHelper.extractTumorType(cpctEcrfModel, sample);
        final Double tumorPercentage = limsModel.tumorPercentageForSample(sample);
        final String sampleRecipient = centerModel.getAddresseeStringForSample(sample);
        final SampleReport sampleReport = ImmutableSampleReport.of(sample, tumorType, tumorPercentage, limsModel.barcodeForSample(sample),
                limsModel.bloodBarcodeForSample(sample), limsModel.arrivalDateForSample(sample),
                limsModel.bloodArrivalDateForSample(sample), limsModel.labProceduresForSample(sample), sampleRecipient);
        final NotSequencedPatientReport report =
                ImmutableNotSequencedPatientReport.of(sampleReport, reason, study, Optional.ofNullable(comments), signaturePath);
        reportWriter.writeNonSequenceableReport(report);
    }
}
