package com.hartwig.hmftools.patientreporter.algo;

import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientreporter.BaseReporterData;
import com.hartwig.hmftools.patientreporter.ImmutableNotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.NotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class NotAnalysableReporter {

    @NotNull
    abstract BaseReporterData baseReporterData();

    public NotAnalysedPatientReport run(@NotNull final String sample, @NotNull final NotAnalysableReason reason,
            @Nullable final String comments) {
        final NotAnalysableStudy study = NotAnalysableStudy.fromSample(sample);
        assert study != null;

        final PatientTumorLocation patientTumorLocation =
                PatientReporterHelper.extractPatientTumorLocation(baseReporterData().patientTumorLocations(), sample);
        final Lims lims = baseReporterData().limsModel();
        final Double tumorPercentage = lims.tumorPercentageForSample(sample);
        final String sampleRecipient = baseReporterData().centerModel().getAddresseeStringForSample(sample);

        final SampleReport sampleReport = ImmutableSampleReport.of(sample,
                patientTumorLocation,
                tumorPercentage,
                lims.arrivalDateForSample(sample),
                null,
                lims.labProceduresForSample(sample),
                sampleRecipient);

        return ImmutableNotAnalysedPatientReport.of(sampleReport,
                reason,
                study,
                Optional.ofNullable(comments),
                baseReporterData().signaturePath());
    }
}
