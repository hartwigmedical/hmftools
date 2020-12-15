package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;

import org.jetbrains.annotations.NotNull;

public class PatientReportUtils {

    @NotNull
    public static LimsCohortConfigData buildTestCohortModel(@NotNull String cohortId, boolean hospitalCentraId, boolean reportGermline,
            boolean reportGermlineFlag, boolean reportConclusion, boolean reportViral, boolean requireHospitalId,
            boolean requireHospitalPAId, boolean requireHospitalPersonsStudy, boolean requireHospitalPersonsRequester,
            boolean requirePatientIdForPdfName, boolean requireSubmissionInformation, boolean requireAdditionalInfromationForSidePanel) {
        return ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortId)
                .hospitalCentraId(hospitalCentraId)
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
                .requireAdditionalInfromationForSidePanel(requireAdditionalInfromationForSidePanel)
                .build();
    }
}
