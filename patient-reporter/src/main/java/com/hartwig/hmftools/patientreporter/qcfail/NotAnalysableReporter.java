package com.hartwig.hmftools.patientreporter.qcfail;

import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientreporter.BaseReportData;
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
    abstract BaseReportData baseReportData();

    public NotAnalysedPatientReport run(@NotNull final String sample, @NotNull final NotAnalysableReason reason,
            @Nullable final String comments) {
        final Lims lims = baseReportData().limsModel();

        final NotAnalysableStudy study = NotAnalysableStudy.fromSample(sample);

        assert study != null;

        final PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(baseReportData().patientTumorLocations(), sample);

        boolean isCoreSample = baseReportData().limsModel().isCoreSample(sample);
        final SampleReport sampleReport = ImmutableSampleReport.of(sample,
                patientTumorLocation,
                lims.purityShallowSeq(sample, false),
                lims.arrivalDateForSample(sample),
                null,
                lims.labProceduresForSample(sample),
                isCoreSample
                        ? baseReportData().centerModel().getCoreRecipients(lims.projectNameDVO(sample))
                        : baseReportData().centerModel().getAddresseeStringForSample(sample),
                lims.labelSample(sample),
                lims.projectNameDVO(sample),
                isCoreSample ? lims.contactEmail(sample) : "",
                isCoreSample ? lims.contactName(sample) : "",
                lims.patientNumber(sample),
                lims.isCoreSample(sample));

        return ImmutableNotAnalysedPatientReport.of(sampleReport,
                reason,
                study,
                Optional.ofNullable(comments),
                baseReportData().signaturePath(),
                baseReportData().logoRVAPath());
    }
}
