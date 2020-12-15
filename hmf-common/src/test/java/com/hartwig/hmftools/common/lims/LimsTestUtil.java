package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;

import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LimsTestUtil {

    private LimsTestUtil() {
    }

    @NotNull
    static LocalDate toDate(@NotNull String date) {
        return LocalDate.parse(date, LimsConstants.DATE_FORMATTER);
    }

    @NotNull
    public static ImmutableLimsJsonSampleData.Builder createLimsSampleDataBuilder() {
        return ImmutableLimsJsonSampleData.builder()
                .sampleId(Strings.EMPTY)
                .patientId(Strings.EMPTY)
                .tumorBarcode(Strings.EMPTY)
                .refBarcode(Strings.EMPTY)
                .arrivalDate(Strings.EMPTY)
                .dnaConcentration(Strings.EMPTY)
                .primaryTumor(Strings.EMPTY)
                .labSopVersions(Strings.EMPTY)
                .submission(Strings.EMPTY)
                .germlineReportingLevel(Strings.EMPTY)
                .reportGermlineVariants(false)
                .shallowSeq(false)
                .reportViralInsertions(false)
                .cohort(Strings.EMPTY)
                .analysisTypeSample(Strings.EMPTY);
    }

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
