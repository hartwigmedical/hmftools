package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;

import org.jetbrains.annotations.NotNull;

public final class PatientReportUtils {

    private PatientReportUtils() {
    }

    @NotNull
    public static LimsCohortConfig buildTestCohortModel(@NotNull String cohortId, boolean hospitalCenterId, boolean reportGermline,
            boolean reportGermlineFlag, boolean reportConclusion, boolean reportViral, boolean requireHospitalId,
            boolean requireHospitalPAId, boolean requireHospitalPersonsStudy, boolean requireHospitalPersonsRequester,
            boolean requirePatientIdForPdfName, boolean requireSubmissionInformation, boolean requireAdditionalInformationForSidePanel) {
        return ImmutableLimsCohortConfig.builder()
                .cohortId(cohortId)
                .hospitalCenterId(hospitalCenterId)
                .reportGermline(reportGermline)
                .reportGermlineFlag(reportGermlineFlag)
                .reportConclusion(reportConclusion)
                .reportViral(reportViral)
                .requireHospitalId(requireHospitalId)
                .requireHospitalPAId(requireHospitalPAId)
                .requireHospitalPersonsStudy(requireHospitalPersonsStudy)
                .requireHospitalPersonsRequester(requireHospitalPersonsRequester)
                .requirePatientIdForPdfName(requirePatientIdForPdfName)
                .requireSubmissionInformation(requireSubmissionInformation)
                .requireAdditionalInformationForSidePanel(requireAdditionalInformationForSidePanel)
                .build();
    }
}
