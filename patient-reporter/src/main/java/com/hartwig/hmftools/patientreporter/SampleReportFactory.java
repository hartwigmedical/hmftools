package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsCohort;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SampleReportFactory {

    private static final Logger LOGGER = LogManager.getLogger(SampleReportFactory.class);

    private SampleReportFactory() {
    }

    @NotNull
    public static SampleReport fromLimsModel(@NotNull SampleMetadata sampleMetadata, @NotNull Lims lims,
            @Nullable PatientPrimaryTumor patientPrimaryTumor) {
        String refSampleBarcode = sampleMetadata.refSampleBarcode();
        String refSampleId = sampleMetadata.refSampleId();
        String tumorSampleBarcode = sampleMetadata.tumorSampleBarcode();
        String tumorSampleId = sampleMetadata.tumorSampleId();

        LocalDate arrivalDateRefSample = null;

        if (refSampleBarcode != null && refSampleId != null) {
            lims.validateSampleBarcodeCombination(refSampleBarcode, refSampleId, tumorSampleBarcode, tumorSampleId);

            arrivalDateRefSample = lims.arrivalDate(refSampleBarcode, refSampleId);
            if (arrivalDateRefSample == null) {
                LOGGER.warn("Could not find arrival date for ref sample: {}", refSampleId);
            }
        }

        LocalDate arrivalDateTumorSample = lims.arrivalDate(tumorSampleBarcode, tumorSampleId);
        if (arrivalDateTumorSample == null) {
            LOGGER.warn("Could not find arrival date for tumor sample: {}", tumorSampleId);
        }

        String hospitalPathologySampleId = lims.hospitalPathologySampleId(tumorSampleBarcode);

        LimsCohortConfigData cohortdata = lims.cohortConfig(tumorSampleBarcode);

        LOGGER.info("Cohort type of this sample is: {}", cohortdata.cohortId());

        String hospitalPatientId =
                checkHospitalPatientId(lims.hospitalPatientId(tumorSampleBarcode), tumorSampleId, cohortdata);

        return ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientPrimaryTumor(patientPrimaryTumor)
                .germlineReportingLevel(lims.germlineReportingChoice(tumorSampleBarcode))
                .reportViralInsertions(lims.reportViralInsertions(tumorSampleBarcode))
                .refArrivalDate(arrivalDateRefSample)
                .tumorArrivalDate(arrivalDateTumorSample)
                .shallowSeqPurityString(lims.purityShallowSeq(tumorSampleBarcode))
                .labProcedures(lims.labProcedures(tumorSampleBarcode))
                .cohort(lims.cohort(tumorSampleBarcode))
                .projectName(lims.projectName(tumorSampleBarcode))
                .submissionId(lims.submissionId(tumorSampleBarcode))
                .hospitalContactData(lims.hospitalContactData(tumorSampleBarcode))
                .hospitalPatientId(hospitalPatientId)
                .hospitalPathologySampleId(toHospitalPathologySampleIdForReport(hospitalPathologySampleId, tumorSampleId, cohortdata))
                .build();
    }

    @VisibleForTesting
    static String checkHospitalPatientId(@NotNull String hospitalPatientId, @NotNull String sampleId,
            @NotNull LimsCohortConfigData cohortdata) {
        if (cohortdata.requireHospitalId()) {
            if (hospitalPatientId.equals(Lims.NOT_AVAILABLE_STRING) || hospitalPatientId.equals(Strings.EMPTY)) {
                LOGGER.warn("Missing hospital patient sample ID for sample '{}': {}. Please fix!", sampleId, hospitalPatientId);
            }
        }
        return hospitalPatientId;
    }

    @VisibleForTesting
    @Nullable
    static String toHospitalPathologySampleIdForReport(@NotNull String hospitalPathologySampleId, @NotNull String tumorSampleId,
            @Nullable LimsCohortConfigData cohortdata) {
        if (cohortdata.requireHospitalPAId()) {
            if (!hospitalPathologySampleId.equals(Lims.NOT_AVAILABLE_STRING) && !hospitalPathologySampleId.isEmpty()
                    && isValidHospitalPathologySampleId(hospitalPathologySampleId)) {
                return hospitalPathologySampleId;
            } else {

                LOGGER.warn("Missing or invalid hospital pathology sample ID for sample '{}': {}. Please fix!",
                        tumorSampleId,
                        hospitalPathologySampleId);

                return null;
            }
        } else {
            if (!hospitalPathologySampleId.isEmpty() && !hospitalPathologySampleId.equals(Lims.NOT_AVAILABLE_STRING)) {
                LOGGER.info("Skipping hospital pathology sample ID for sample '{}': {}", hospitalPathologySampleId, tumorSampleId);
            }

            return null;
        }
    }

    private static boolean isValidHospitalPathologySampleId(@NotNull String hospitalPathologySampleId) {
        boolean tMatch = hospitalPathologySampleId.startsWith("T") && hospitalPathologySampleId.substring(1, 3).matches("[0-9]+")
                && hospitalPathologySampleId.substring(3, 4).equals("-") && hospitalPathologySampleId.substring(4, 9).matches("[0-9]+");

        boolean cMatch = hospitalPathologySampleId.startsWith("C") && hospitalPathologySampleId.substring(1, 3).matches("[0-9]+")
                && hospitalPathologySampleId.substring(3, 4).equals("-") && hospitalPathologySampleId.substring(4, 9).matches("[0-9]+");

        return tMatch || cMatch;
    }
}
