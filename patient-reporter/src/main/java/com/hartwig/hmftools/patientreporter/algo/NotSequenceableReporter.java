package com.hartwig.hmftools.patientreporter.algo;

import java.io.FileNotFoundException;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.patientreporter.lims.TumorPercentages;
import com.hartwig.hmftools.patientreporter.report.ReportWriter;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.exception.DRException;

public class NotSequenceableReporter {

    @NotNull
    private final CpctEcrfModel cpctEcrfModel;
    @NotNull
    private final ReportWriter reportWriter;
    @NotNull
    private final TumorPercentages tumorPercentages;

    public NotSequenceableReporter(@NotNull final CpctEcrfModel cpctEcrfModel,
            @NotNull final ReportWriter reportWriter,
            @NotNull final TumorPercentages tumorPercentages) {
        this.cpctEcrfModel = cpctEcrfModel;
        this.reportWriter = reportWriter;
        this.tumorPercentages = tumorPercentages;
    }

    public void run(@NotNull final String sample, @NotNull final NotSequenceableReason reason)
            throws FileNotFoundException, DRException {
        final String tumorType = PatientReporterHelper.extractTumorType(cpctEcrfModel, sample);
        final String tumorPercentageString =
                Long.toString(Math.round(100D * tumorPercentages.findTumorPercentageForSample(sample))) + "%";
        reportWriter.writeNonSequenceableReport(sample, tumorType, tumorPercentageString, reason);
    }
}
