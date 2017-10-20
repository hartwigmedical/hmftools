package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.util.Optional;

import com.hartwig.hmftools.common.lims.LimsJsonModel;
import com.hartwig.hmftools.patientreporter.BaseReporterData;
import com.hartwig.hmftools.patientreporter.ImmutableNotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.NotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.exception.DRException;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class NotSequenceableReporter {

    @NotNull
    abstract BaseReporterData baseReporterData();

    public NotSequencedPatientReport run(@NotNull final String sample, @NotNull final NotSequenceableReason reason,
            @Nullable final String comments) throws IOException, DRException {
        final NotSequenceableStudy study = NotSequenceableStudy.fromSample(sample);
        assert study != null;
        final String tumorType = PatientReporterHelper.extractTumorType(baseReporterData().cpctEcrfModel(), sample);
        final LimsJsonModel limsModel = baseReporterData().limsModel();
        final Double tumorPercentage = limsModel.tumorPercentageForSample(sample);
        final String sampleRecipient = baseReporterData().centerModel().getAddresseeStringForSample(sample);
        final SampleReport sampleReport =
                ImmutableSampleReport.of(sample, tumorType, tumorPercentage, limsModel.arrivalDateForSample(sample),
                        limsModel.bloodArrivalDateForSample(sample), limsModel.labProceduresForSample(sample), sampleRecipient);
        return ImmutableNotSequencedPatientReport.of(sampleReport, reason, study, Optional.ofNullable(comments),
                baseReporterData().signaturePath());
    }
}
