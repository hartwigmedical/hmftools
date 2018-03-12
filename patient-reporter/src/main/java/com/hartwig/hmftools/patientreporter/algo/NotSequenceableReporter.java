package com.hartwig.hmftools.patientreporter.algo;

import java.util.Optional;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientreporter.BaseReporterData;
import com.hartwig.hmftools.patientreporter.ImmutableNotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.NotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class NotSequenceableReporter {

    @NotNull
    abstract BaseReporterData baseReporterData();

    public NotSequencedPatientReport run(@NotNull final String sample, @NotNull final NotSequenceableReason reason,
            @Nullable final String comments) {
        final NotSequenceableStudy study = NotSequenceableStudy.fromSample(sample);
        assert study != null;
        final String cancerType = PatientReporterHelper.extractCancerType(baseReporterData().patientsCancerTypes(), sample);
        final Lims lims = baseReporterData().limsModel();
        final Double tumorPercentage = lims.tumorPercentageForSample(sample);
        final String sampleRecipient = baseReporterData().centerModel().getAddresseeStringForSample(sample);
        final SampleReport sampleReport = ImmutableSampleReport.of(sample, cancerType,
                tumorPercentage,
                lims.arrivalDateForSample(sample),
                null,
                lims.labProceduresForSample(sample),
                sampleRecipient);
        return ImmutableNotSequencedPatientReport.of(sampleReport,
                reason,
                study,
                Optional.ofNullable(comments),
                baseReporterData().signaturePath());
    }
}
