package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;
import java.time.format.DateTimeParseException;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.HospitalModel;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Lims {

    private static final Logger LOGGER = LogManager.getLogger(Lims.class);

    public static final String NOT_AVAILABLE_STRING = "N/A";
    public static final String NOT_PERFORMED_STRING = "not performed";
    public static final String PURITY_NOT_RELIABLE_STRING = "below detection threshold";

    @NotNull
    private final LimsCohortModel limsCohortModel;
    @NotNull
    private final HospitalModel hospitalModel;
    @NotNull
    private final Map<String, LimsJsonSampleData> dataPerSampleBarcode;
    @NotNull
    private final Map<String, LimsJsonSubmissionData> dataPerSubmission;
    @NotNull
    private final Map<String, LimsShallowSeqData> shallowSeqPerSampleBarcode;
    @NotNull
    private final Map<String, LocalDate> preLimsArrivalDatesPerSampleId;
    @NotNull
    private final Set<String> samplesWithoutSamplingDate;
    @NotNull
    private final Set<String> patientsWithoutCuratedPrimaryTumor;
    @NotNull
    private final Set<String> blacklistedPatients;

    public Lims(@NotNull final LimsCohortModel limsCohortModel, @NotNull final HospitalModel hospitalModel,
            @NotNull final Map<String, LimsJsonSampleData> dataPerSampleBarcode,
            @NotNull final Map<String, LimsJsonSubmissionData> dataPerSubmission,
            @NotNull final Map<String, LimsShallowSeqData> shallowSeqPerSampleBarcode,
            @NotNull final Map<String, LocalDate> preLimsArrivalDatesPerSampleId, @NotNull final Set<String> samplesWithoutSamplingDate,
            @NotNull final Set<String> patientsWithoutCuratedPrimaryTumor, @NotNull final Set<String> blacklistedPatients) {
        this.limsCohortModel = limsCohortModel;
        this.hospitalModel = hospitalModel;
        this.dataPerSampleBarcode = dataPerSampleBarcode;
        this.dataPerSubmission = dataPerSubmission;
        this.shallowSeqPerSampleBarcode = shallowSeqPerSampleBarcode;
        this.preLimsArrivalDatesPerSampleId = preLimsArrivalDatesPerSampleId;
        this.samplesWithoutSamplingDate = samplesWithoutSamplingDate;
        this.patientsWithoutCuratedPrimaryTumor = patientsWithoutCuratedPrimaryTumor;
        this.blacklistedPatients = blacklistedPatients;
    }

    public int sampleBarcodeCount() {
        return dataPerSampleBarcode.size();
    }

    @NotNull
    public Set<String> sampleBarcodes() {
        return dataPerSampleBarcode.keySet();
    }

    public void validateSampleBarcodeCombination(@NotNull String refBarcode, @NotNull String refSampleId, @NotNull String tumorBarcode,
            @NotNull String tumorSampleId) {
        LimsJsonSampleData tumorSampleData = dataPerSampleBarcode.get(tumorBarcode);
        if (tumorSampleData == null) {
            LOGGER.warn("Could not find entry for barcode '{}' in LIMS.", tumorBarcode);
        } else {
            if (!tumorSampleData.sampleId().equals(tumorSampleId)) {
                LOGGER.warn("Mismatching tumor sample name. Provided='{}'. LIMS='{}'", tumorSampleId, tumorSampleData.sampleId());
            }

            String limsRefBarcode = tumorSampleData.refBarcode();
            if (limsRefBarcode == null || !limsRefBarcode.equals(refBarcode)) {
                LOGGER.warn("Mismatching ref sample barcode. Provided='{}'. LIMS='{}'", refBarcode, limsRefBarcode);
            } else {
                LimsJsonSampleData refSampleData = dataPerSampleBarcode.get(limsRefBarcode);
                if (refSampleData == null) {
                    LOGGER.warn("No ref sample data for sample with barcode '{}'", limsRefBarcode);
                } else if (!refSampleData.sampleId().equals(refSampleId)) {
                    LOGGER.warn("Mismatching ref sample name. Provided='{}'. LIMS='{}'", refSampleId, refSampleData.sampleId());
                }
            }
        }
    }

    @NotNull
    public String patientId(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        return sampleData != null ? sampleData.patientId() : NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String sampleId(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        return sampleData != null ? sampleData.sampleId() : NOT_AVAILABLE_STRING;
    }

    @Nullable
    public LocalDate arrivalDate(@NotNull String sampleBarcode, @NotNull String fallbackSampleId) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        LocalDate arrivalDate = sampleData != null ? getNullableDate(sampleData.arrivalDate()) : null;

        if (arrivalDate == null) {
            arrivalDate = preLimsArrivalDatesPerSampleId.get(fallbackSampleId);
        }

        return arrivalDate;
    }

    @Nullable
    public LocalDate samplingDate(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            String samplingDateString = sampleData.samplingDate();
            return getNullableDate(samplingDateString);
        }
        return null;
    }

    public boolean confirmedToHaveNoSamplingDate(@NotNull String sampleId) {
        return samplesWithoutSamplingDate.contains(sampleId);
    }

    public boolean isBlacklisted(@NotNull String patientId) {
        return blacklistedPatients.contains(patientId);
    }

    public boolean needsNoCuratedPrimaryTumor(@NotNull String patientId) {
        return patientsWithoutCuratedPrimaryTumor.contains(patientId);
    }

    @NotNull
    public String submissionId(@NotNull String sampleBarcode) {
        String submission = submission(sampleBarcode);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.submission() : NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String projectName(@NotNull String sampleBarcode) {
        String submission = submission(sampleBarcode);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.projectName() : NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String requesterName(@NotNull String sampleBarcode) {
        String submission = submission(sampleBarcode);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.reportContactName() : NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String requesterEmail(@NotNull String sampleBarcode) {
        String submission = submission(sampleBarcode);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.reportContactEmail() : NOT_AVAILABLE_STRING;
    }

    @Nullable
    public Integer dnaNanograms(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            try {
                // LIMS stores the amount of nanograms per micro liter.
                return (int) Math.round(Double.parseDouble(sampleData.dnaConcentration()) * LimsConstants.DNA_MICRO_LITERS);
            } catch (NumberFormatException e) {
                return null;
            }
        }
        return null;
    }

    @NotNull
    public String purityShallowSeq(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            boolean purityShallowExecuted = shallowSeqExecuted(sampleBarcode);
            LimsShallowSeqData shallowSeq = shallowSeqPerSampleBarcode.get(sampleBarcode);

            if (purityShallowExecuted && shallowSeq == null) {
                LOGGER.warn("BFX lims and lab status on shallow seq do not match for sample '{}'!", sampleBarcode);
            } else {
                if (purityShallowExecuted) {
                    if (!shallowSeq.hasReliablePurity()) {
                        return PURITY_NOT_RELIABLE_STRING;
                    } else {
                        try {
                            return Math.round(Double.parseDouble(shallowSeq.purityShallowSeq()) * 100) + "%";
                        } catch (NumberFormatException e) {
                            LOGGER.warn("Could not convert shallow seq to a percentage '{}'", shallowSeq.purityShallowSeq());
                            return NOT_AVAILABLE_STRING;
                        }
                    }
                } else {
                    return NOT_PERFORMED_STRING;
                }
            }
        }
        return NOT_PERFORMED_STRING;
    }

    @NotNull
    public String pathologyTumorPercentage(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            String tumorPercentageString = sampleData.pathologyTumorPercentage();
            if (tumorPercentageString == null) {
                return NOT_AVAILABLE_STRING;
            } else {
                try {
                    return Math.round(Double.parseDouble(tumorPercentageString)) + "%";
                } catch (NumberFormatException e) {
                    return NOT_AVAILABLE_STRING;
                }
            }
        }
        return NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String primaryTumor(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            return sampleData.primaryTumor();
        }
        // No warning raised since initially this information was not tracked so this will be missing for early samples.
        return NOT_AVAILABLE_STRING;
    }

    @Nullable
    public String biopsyLocation(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            return sampleData.biopsySite();
        }
        // No warning raised since initially this information was not tracked so this will be missing for early samples.
        return null;
    }

    @NotNull
    public String labProcedures(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            return sampleData.labProcedures();
        }
        LOGGER.warn("Could not find lab procedures for sample '{}' in LIMS", sampleBarcode);
        return NOT_AVAILABLE_STRING;
    }

    @Nullable
    public LimsCohortConfig cohortConfig(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        String cohortString = sampleData != null ? sampleData.cohort() : null;

        return limsCohortModel.queryCohortData(cohortString, sampleId(sampleBarcode));
    }

    @Nullable
    public LimsAnalysisType extractAnalysisType(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            String analysisTypeString = sampleData.analysisType();
            return analysisTypeString != null ? LimsAnalysisType.extractAnalysisType(analysisTypeString) : null;
        } else {
            return null;
        }
    }

    @NotNull
    public HospitalContactData hospitalContactData(@NotNull String sampleBarcode) {
        String sampleId = sampleId(sampleBarcode);
        HospitalContactData data = hospitalModel.queryHospitalData(sampleId,
                cohortConfig(sampleBarcode),
                requesterName(sampleBarcode),
                requesterEmail(sampleBarcode));

        if (data == null) {
            LOGGER.warn("Could not find hospital data for sample '{}' with barcode '{}'", sampleId, sampleBarcode);
            data = ImmutableHospitalContactData.builder()
                    .hospitalPI(NOT_AVAILABLE_STRING)
                    .requesterName(NOT_AVAILABLE_STRING)
                    .requesterEmail(NOT_AVAILABLE_STRING)
                    .hospitalName(NOT_AVAILABLE_STRING)
                    .hospitalAddress(NOT_AVAILABLE_STRING)
                    .build();
        }

        return data;
    }

    @NotNull
    public String hospitalPatientId(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            String hospitalPatientId = sampleData.hospitalPatientId();
            return hospitalPatientId != null ? hospitalPatientId : NOT_AVAILABLE_STRING;
        }
        return NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String hospitalPathologySampleId(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            String hospitalPathologySampleId = sampleData.hospitalPathologySampleId();
            return hospitalPathologySampleId != null ? hospitalPathologySampleId : NOT_AVAILABLE_STRING;
        }
        return NOT_AVAILABLE_STRING;
    }

    @NotNull
    public LimsGermlineReportingLevel germlineReportingChoice(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        if (sampleData != null) {
            String germlineReportingLevelString = sampleData.germlineReportingLevel();

            if (cohortConfig(sampleBarcode) != null) {
                if (cohortConfig(sampleBarcode).reportGermline() && germlineReportingLevelString == null) {
                    LOGGER.warn("Germline choice is unknown but is expected");
                }
            }

            return germlineReportingLevelString != null ? LimsGermlineReportingLevel.fromLimsInputs(reportGermlineVariants(sampleBarcode),
                    germlineReportingLevelString,
                    sampleId(sampleBarcode),
                    cohortConfig(sampleBarcode)) : LimsGermlineReportingLevel.NO_REPORTING;
        } else {
            return LimsGermlineReportingLevel.NO_REPORTING;
        }
    }

    @VisibleForTesting
    boolean reportGermlineVariants(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        LimsCohortConfig cohort = cohortConfig(sampleBarcode);

        return LimsChecker.checkGermlineVariants(sampleData, cohort, sampleId(sampleBarcode));
    }

    public boolean reportViralPresence(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        LimsCohortConfig cohort = cohortConfig(sampleBarcode);

        return LimsChecker.checkViralInsertions(sampleData, cohort, sampleId(sampleBarcode));
    }

    public boolean reportPgx(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        LimsCohortConfig cohort = cohortConfig(sampleBarcode);

        return LimsChecker.checkReportPgx(sampleData, cohort, sampleId(sampleBarcode));
    }

    @Nullable
    private String submission(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        return sampleData != null ? sampleData.submission() : null;
    }

    private boolean shallowSeqExecuted(@NotNull String sampleBarcode) {
        LimsJsonSampleData sampleData = dataPerSampleBarcode.get(sampleBarcode);
        assert sampleData != null;

        return sampleData.shallowSeq();
    }

    @Nullable
    private static LocalDate getNullableDate(@Nullable String dateString) {
        if (dateString == null) {
            return null;
        }

        try {
            return LocalDate.parse(dateString, LimsConstants.DATE_FORMATTER);
        } catch (DateTimeParseException e) {
            return null;
        }
    }
}
