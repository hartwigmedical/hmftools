package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsCoreCohort;
import com.hartwig.hmftools.common.lims.LimsStudy;

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
            @Nullable PatientTumorLocation patientTumorLocation) {
        String refSampleBarcode = sampleMetadata.refSampleBarcode();
        String refSampleId = sampleMetadata.refSampleId();
        String tumorSampleBarcode = sampleMetadata.tumorSampleBarcode();
        String tumorSampleId = sampleMetadata.tumorSampleId();

        lims.validateSampleBarcodeCombination(refSampleBarcode, refSampleId, tumorSampleBarcode, tumorSampleId);

        LocalDate arrivalDateRefSample = lims.arrivalDate(refSampleBarcode, refSampleId);
        if (arrivalDateRefSample == null) {
            LOGGER.warn("Could not find arrival date for ref sample: {}", refSampleId);
        }

        LocalDate arrivalDateTumorSample = lims.arrivalDate(tumorSampleBarcode, tumorSampleId);
        if (arrivalDateTumorSample == null) {
            LOGGER.warn("Could not find arrival date for tumor sample: {}", tumorSampleId);
        }

        String hospitalPathologySampleId = lims.hospitalPathologySampleId(tumorSampleBarcode);

        String cohort = lims.cohort(tumorSampleBarcode);
        LimsStudy type = LimsStudy.fromSampleId(tumorSampleId);
        LimsCoreCohort coreCohort = LimsCoreCohort.fromSampleId(tumorSampleId);

        if (cohort.isEmpty()) {
            if (coreCohort.equals(LimsCoreCohort.NON_CORE)) {
                cohort = type.toString();
            }
        }

        return ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientTumorLocation(patientTumorLocation)
                .refArrivalDate(arrivalDateRefSample)
                .tumorArrivalDate(arrivalDateTumorSample)
                .shallowSeqPurityString(lims.purityShallowSeq(tumorSampleBarcode))
                .labProcedures(lims.labProcedures(tumorSampleBarcode))
                .cohort(cohort)
                .projectName(lims.projectName(tumorSampleBarcode))
                .submissionId(lims.submissionId(tumorSampleBarcode))
                .hospitalContactData(lims.hospitalContactData(tumorSampleBarcode))
                .hospitalPatientId(lims.hospitalPatientId(tumorSampleBarcode))
                .hospitalPathologySampleId(toHospitalPathologySampleIdForReport(hospitalPathologySampleId, tumorSampleId))
                .build();
    }

    @VisibleForTesting
    @Nullable
    static String toHospitalPathologySampleIdForReport(@NotNull String hospitalPathologySampleId, @NotNull String tumorSampleId) {
        LimsStudy study = LimsStudy.fromSampleId(tumorSampleId);

        if (study == LimsStudy.CORE || study == LimsStudy.WIDE) {
            if (!hospitalPathologySampleId.equals(Lims.NOT_AVAILABLE_STRING) && !hospitalPathologySampleId.isEmpty()
                    && isValidHospitalPathologySampleId(hospitalPathologySampleId)) {
                return hospitalPathologySampleId;
            } else {
                if (study == LimsStudy.WIDE) {
                    LOGGER.warn("Missing or invalid hospital pathology sample ID for sample '{}': {}. Please fix!",
                            tumorSampleId,
                            hospitalPathologySampleId);
                } else {
                    if (!hospitalPathologySampleId.isEmpty()) {
                        LOGGER.warn("No valid hospital pathology sample ID found for '{}': {}", tumorSampleId, hospitalPathologySampleId);
                    }
                }
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
