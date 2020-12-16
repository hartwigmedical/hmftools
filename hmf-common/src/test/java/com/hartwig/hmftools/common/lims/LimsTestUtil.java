package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;

import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;

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
    public static LimsCohortConfig createAllDisabledCohortConfig(@NotNull String cohortId) {
        return allDisabledBuilder().cohortId(cohortId).build();
    }

    @NotNull
    public static LimsCohortConfig createConfigForHospitalModel(@NotNull String cohortId, boolean requireHospitalPersonsStudy,
            boolean requireHospitalPersonsRequester) {
        return allDisabledBuilder().cohortId(cohortId)
                .sampleContainsHospitalCenterId(true)
                .requireHospitalPersonsStudy(requireHospitalPersonsStudy)
                .requireHospitalPersonsRequester(requireHospitalPersonsRequester)
                .build();
    }

    @NotNull
    public static LimsCohortConfig createTestCohortConfig(@NotNull String cohortId, boolean sampleContainsHospitalCenterId,
            boolean reportGermline, boolean reportGermlineFlag, boolean reportConclusion, boolean reportViral, boolean requireHospitalId,
            boolean requireHospitalPAId, boolean requireHospitalPersonsStudy, boolean requireHospitalPersonsRequester,
            boolean requirePatientIdForPdfName, boolean requireSubmissionInformation, boolean requireAdditionalInformationForSidePanel) {
        return ImmutableLimsCohortConfig.builder()
                .cohortId(cohortId)
                .sampleContainsHospitalCenterId(sampleContainsHospitalCenterId)
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

    @NotNull
    private static ImmutableLimsCohortConfig.Builder allDisabledBuilder() {
        return ImmutableLimsCohortConfig.builder()
                .sampleContainsHospitalCenterId(false)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(false)
                .reportViral(false)
                .requireHospitalId(false)
                .requireHospitalPAId(false)
                .requireHospitalPersonsStudy(false)
                .requireHospitalPersonsRequester(false)
                .requirePatientIdForPdfName(false)
                .requireSubmissionInformation(false)
                .requireAdditionalInformationForSidePanel(false);
    }
}
