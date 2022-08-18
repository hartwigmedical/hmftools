package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsChecker;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SampleReportFactory {

    private static final Logger LOGGER = LogManager.getLogger(SampleReportFactory.class);

    private SampleReportFactory() {
    }

    @Nullable
    public static String interpretRefBarcode(@Nullable String refSampleBarcode){
        String interpretRefSampleBarcode = Strings.EMPTY;
        if (refSampleBarcode != null) {
            if (refSampleBarcode.contains("-")) {
                interpretRefSampleBarcode = refSampleBarcode.split("-")[0];
            } else if (refSampleBarcode.contains("_")) {
                interpretRefSampleBarcode = refSampleBarcode.split("_")[0];
            } else {
                interpretRefSampleBarcode = refSampleBarcode;
            }
        } else {
            interpretRefSampleBarcode = refSampleBarcode;
        }
        return interpretRefSampleBarcode;
    }

    @NotNull
    public static SampleReport fromLimsModel(@NotNull SampleMetadata sampleMetadata, @NotNull Lims lims,
            @Nullable PatientPrimaryTumor patientPrimaryTumor, boolean allowDefaultCohortConfig) {
        String interpretRefSampleBarcode = interpretRefBarcode(sampleMetadata.refSampleBarcode());
        String refSampleId = sampleMetadata.refSampleId();
        String tumorSampleBarcode = sampleMetadata.tumorSampleBarcode();
        String tumorSampleId = sampleMetadata.tumorSampleId();

        LocalDate arrivalDateRefSample = null;

        if (interpretRefSampleBarcode != null && refSampleId != null) {
            if (!interpretRefSampleBarcode.equals(refSampleId) || !tumorSampleBarcode.equals(tumorSampleId)){
                // Don't need to check for anonymised runs
                lims.validateSampleBarcodeCombination(interpretRefSampleBarcode, refSampleId, tumorSampleBarcode, tumorSampleId);
            }

            arrivalDateRefSample = lims.arrivalDate(interpretRefSampleBarcode, refSampleId);
            if (arrivalDateRefSample == null) {
                LOGGER.warn("Could not find arrival date for ref sample: {}", refSampleId);
            }
        }

        LocalDate arrivalDateTumorSample = lims.arrivalDate(tumorSampleBarcode, tumorSampleId);
        if (arrivalDateTumorSample == null) {
            LOGGER.warn("Could not find arrival date for tumor sample: {}", tumorSampleId);
        }

        String hospitalPathologySampleId = lims.hospitalPathologySampleId(tumorSampleBarcode);

        LimsCohortConfig cohortConfig = lims.cohortConfig(tumorSampleBarcode);
        if (cohortConfig == null) {
            if (allowDefaultCohortConfig) {
                LOGGER.warn("Using DEFAULT cohort config for tumor sample (non-production only!): {}", tumorSampleId);
                cohortConfig = buildDefaultCohortConfig();
            } else {
                throw new IllegalStateException(
                        "Cohort not configured in LIMS for sample '" + tumorSampleId + "' with barcode " + tumorSampleBarcode);
            }
        }

        String hospitalPatientId = lims.hospitalPatientId(tumorSampleBarcode);
        LimsChecker.checkHospitalPatientId(hospitalPatientId, tumorSampleId, cohortConfig);
        String biopsyLocation = lims.biopsyLocation(tumorSampleBarcode);
        String curatedBiopsyLocation = curateBiopsyLocation(biopsyLocation);

        return ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientPrimaryTumor(patientPrimaryTumor)
                .biopsyLocation(curatedBiopsyLocation)
                .germlineReportingLevel(lims.germlineReportingChoice(tumorSampleBarcode, allowDefaultCohortConfig))
                .reportViralPresence(allowDefaultCohortConfig || lims.reportViralPresence(tumorSampleBarcode))
                .reportPharmogenetics(allowDefaultCohortConfig || lims.reportPgx(tumorSampleBarcode))
                .refArrivalDate(arrivalDateRefSample)
                .tumorArrivalDate(arrivalDateTumorSample)
                .shallowSeqPurityString(lims.purityShallowSeq(tumorSampleBarcode))
                .labProcedures(lims.labProcedures(tumorSampleBarcode))
                .cohort(cohortConfig)
                .projectName(lims.projectName(tumorSampleBarcode))
                .submissionId(lims.submissionId(tumorSampleBarcode))
                .hospitalContactData(lims.hospitalContactData(tumorSampleBarcode))
                .hospitalPatientId(hospitalPatientId)
                .hospitalPathologySampleId(LimsChecker.toHospitalPathologySampleIdForReport(hospitalPathologySampleId,
                        tumorSampleId,
                        cohortConfig))
                .build();
    }

    @Nullable
    public static String curateBiopsyLocation(@Nullable String biopsyLocation) {
        String curated = null;
        if (biopsyLocation != null && biopsyLocation.startsWith("Other (please specify below)")) {
            String[] curatedBiopsyLocation = biopsyLocation.split("_");
            if (curatedBiopsyLocation.length == 2) {
                curated = curatedBiopsyLocation[1];
                curated = curated.substring(0,1).toUpperCase() + curated.substring(1);
            } else if (curatedBiopsyLocation.length == 1) {
                curated = "Other";
            }
        } else {
            if (biopsyLocation != null && !biopsyLocation.equals(Strings.EMPTY)) {
                curated = biopsyLocation;
                curated = curated.substring(0,1).toUpperCase() + curated.substring(1);
            }
        }
        return curated;
    }

    @NotNull
    private static LimsCohortConfig buildDefaultCohortConfig() {
        return ImmutableLimsCohortConfig.builder()
                .cohortId("DEFAULT")
                .sampleContainsHospitalCenterId(false)
                .reportGermline(true)
                .reportGermlineFlag(false)
                .reportConclusion(false)
                .reportViral(true)
                .reportPeach(true)
                .requireHospitalId(false)
                .requireHospitalPAId(false)
                .requireHospitalPersonsStudy(false)
                .requireHospitalPersonsRequester(false)
                .requireAdditionalInformationForSidePanel(false)
                .build();
    }
}
