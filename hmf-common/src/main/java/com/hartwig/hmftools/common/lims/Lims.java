package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;
import java.time.format.DateTimeParseException;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Lims {

    private static final Logger LOGGER = LogManager.getLogger(Lims.class);

    static final String NOT_AVAILABLE_STRING = "N/A";
    static final String NOT_DETERMINED_STRING = "not determined";
    static final String NOT_KNOWN_STRING = "not known";

    @NotNull
    private final Map<String, LimsJsonSampleData> dataPerSample;
    @NotNull
    private final Map<String, LimsJsonSubmissionData> dataPerSubmission;
    @NotNull
    private final Map<String, LocalDate> preLimsArrivalDates;
    @NotNull
    private final Set<String> samplesWithoutSamplingDate;
    @NotNull
    private final Map<String, LimsShallowSeqData> shallowSeqPerSample;

    Lims(@NotNull final Map<String, LimsJsonSampleData> dataPerSample, @NotNull final Map<String, LimsJsonSubmissionData> dataPerSubmission,
            @NotNull final Map<String, LocalDate> preLimsArrivalDates, @NotNull final Set<String> samplesWithoutSamplingDate,
            @NotNull final Map<String, LimsShallowSeqData> shallowSeqPerSample) {
        this.dataPerSample = dataPerSample;
        this.dataPerSubmission = dataPerSubmission;
        this.preLimsArrivalDates = preLimsArrivalDates;
        this.samplesWithoutSamplingDate = samplesWithoutSamplingDate;
        this.shallowSeqPerSample = shallowSeqPerSample;
    }

    public int sampleCount() {
        return dataPerSample.size();
    }

    @NotNull
    public Set<String> sampleIds() {
        return dataPerSample.keySet();
    }

    @NotNull
    public String patientId(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        return sampleData != null ? sampleData.patientId() : NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String tumorBarcode(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        if (sampleData != null) {
            String tumorBarcode = sampleData.tumorBarcode();
            if (tumorBarcode.isEmpty()) {
                return NOT_AVAILABLE_STRING;
            } else {
                return tumorBarcode;
            }
        }
        return NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String refBarcode(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        if (sampleData != null) {
            String refBarcode = sampleData.refBarcode();
            if (refBarcode == null || refBarcode.isEmpty()) {
                return NOT_AVAILABLE_STRING;
            } else {
                return refBarcode;
            }
        }
        return NOT_AVAILABLE_STRING;
    }

    @Nullable
    public LocalDate arrivalDate(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        LocalDate arrivalDate = sampleData != null ? getNullableDate(sampleData.arrivalDate()) : null;

        if (arrivalDate == null) {
            arrivalDate = preLimsArrivalDates.get(sampleId);
        }

        return arrivalDate;
    }

    @Nullable
    public LocalDate samplingDate(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        if (sampleData != null) {
            final String samplingDateString = sampleData.samplingDate();
            return getNullableDate(samplingDateString);
        }
        return null;
    }

    public boolean confirmedToHaveNoSamplingDate(@NotNull String sampleId) {
        return samplesWithoutSamplingDate.contains(sampleId);
    }

    @NotNull
    public String submissionId(@NotNull String sampleId) {
        String submission = submission(sampleId);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.submission() : NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String projectName(@NotNull String sampleId) {
        String submission = submission(sampleId);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.projectName() : NOT_AVAILABLE_STRING;
    }

    @Nullable
    public Integer dnaNanograms(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        if (sampleData != null) {
            try {
                // LIMS stores the amount of nanograms per micro liter.
                return (int) Math.round(Double.parseDouble(sampleData.dnaConcentration()) * LimsConstants.DNA_MICRO_LITERS);
            } catch (final NumberFormatException e) {
                return null;
            }
        }
        return null;
    }

    @NotNull
    public String purityShallowSeq(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);

        if (sampleData != null) {
            boolean purityShallowExecuted = shallowSeqExecuted(sampleId);
            LimsShallowSeqData shallowSeq = shallowSeqPerSample.get(sampleId);

            if (purityShallowExecuted && shallowSeq == null) {
                LOGGER.warn("BFX lims and lab status do not match for sample " + sampleId + "!");
            } else {
                if (purityShallowExecuted) {
                    if (shallowSeq.purityShallowSeq().equals("below detection threshold")) {
                        return "below detection threshold";
                    } else {
                        try {
                            return Math.round(Double.parseDouble(shallowSeq.purityShallowSeq()) * 100) + "%";
                        } catch (final NumberFormatException e) {
                            LOGGER.warn("Could not convert shallow seq to a percentage: " + shallowSeq.purityShallowSeq());
                            return NOT_AVAILABLE_STRING;
                        }
                    }
                } else {
                    return NOT_DETERMINED_STRING;
                }
            }
        }
        return NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String pathologyTumorPercentage(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        if (sampleData != null) {
            // Even if pathology tumor percentage has been determined, we still suppress it in case of shallow seq.
            if (shallowSeqExecuted(sampleId)) {
                return NOT_DETERMINED_STRING;
            }

            String tumorPercentageString = sampleData.pathologyTumorPercentage();
            if (tumorPercentageString == null) {
                return NOT_AVAILABLE_STRING;
            } else {
                try {
                    return Long.toString(Math.round(Double.parseDouble(tumorPercentageString))) + "%";
                } catch (final NumberFormatException e) {
                    return NOT_AVAILABLE_STRING;
                }
            }
        }
        return NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String primaryTumor(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        if (sampleData != null) {
            return sampleData.primaryTumor();
        }
        // No warning raised since initially this information was not tracked so this will be missing for early samples.
        return NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String labProcedures(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        if (sampleData != null) {
            return sampleData.labProcedures();
        }
        LOGGER.warn("Could not find lab SOP versions for sample: " + sampleId + " in LIMS");
        return NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String requesterEmail(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        return sampleData != null ? sampleData.requesterEmail() : NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String requesterName(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        return sampleData != null ? sampleData.requesterName() : NOT_AVAILABLE_STRING;
    }

    @NotNull
    public String hospitalPatientId(@NotNull String sampleId) {
        LimsJsonSampleData sampleData = dataPerSample.get(sampleId);
        if (sampleData != null) {
            String hospitalPatientId = sampleData.hospitalPatientId();
            return hospitalPatientId != null ? hospitalPatientId : NOT_KNOWN_STRING;
        }
        return NOT_KNOWN_STRING;
    }

    @NotNull
    public String hospitalPathologySampleId(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            String hospitalPathologySampleId = sampleData.hospitalPathologySampleId();
            return hospitalPathologySampleId != null ? hospitalPathologySampleId : NOT_KNOWN_STRING;
        }
        return NOT_KNOWN_STRING;
    }

    @NotNull
    public LimsGermlineReportingChoice germlineReportingChoice(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            String germlineReportingChoiceString = sampleData.germlineReportingChoice();
            return germlineReportingChoiceString != null ? LimsGermlineReportingChoice.fromLimsGermlineReportingChoiceString(
                    germlineReportingChoiceString,
                    sample) : LimsGermlineReportingChoice.UNKNOWN;
        } else {
            return LimsGermlineReportingChoice.UNKNOWN;
        }
    }

    @Nullable
    private String submission(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.submission() : null;
    }

    private boolean shallowSeqExecuted(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        assert sampleData != null;

        // TODO: Cleanup once we have switched lims to using "true/false" exclusively
        if (sampleData.shallowSeq().equalsIgnoreCase("1") || sampleData.shallowSeq().equalsIgnoreCase("true")) {
            return true;
        } else {
            String labRemarks = sampleData.labRemarks();
            String labRemarksLowerCase = labRemarks != null ? labRemarks.toLowerCase() : Strings.EMPTY;
            return labRemarksLowerCase.contains("cpctwide") || labRemarksLowerCase.contains("shallowseq");
        }
    }

    @Nullable
    private static LocalDate getNullableDate(@Nullable final String dateString) {
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
