package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;
import java.time.format.DateTimeParseException;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Lims {

    private static final Logger LOGGER = LogManager.getLogger(Lims.class);

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
    public String patientId(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.patientId() : "N/A";
    }

    @NotNull
    public String requesterEmail(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.requesterEmail() : "N/A";
    }

    @NotNull
    public String requesterName(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.requesterName() : "N/A";
    }

    @NotNull
    public String submissionId(@NotNull final String sample) {
        String submission = submission(sample);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.submission() : "N/A";
    }

    @Nullable
    public String hospitalPatientId(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.hospitalPatientId() : null;
    }

    @NotNull
    public String projectName(@NotNull final String sample) {
        String submission = submission(sample);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.projectName() : "N/A";
    }

    @NotNull
    public String barcodeTumor(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            String tumorBarcode = sampleData.tumorBarcodeId();
            if (tumorBarcode.isEmpty()) {
                return "not determined";
            } else {
                return tumorBarcode;
            }
        }
        return "N/A";
    }

    @NotNull
    public String barcodeReference(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            String refBarcode = sampleData.refBarcodeId();
            if (refBarcode == null || refBarcode.isEmpty()) {
                return "not determined";
            } else {
                return refBarcode;
            }
        }
        return "N/A";
    }

    @Nullable
    public LocalDate arrivalDate(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        LocalDate arrivalDate = sampleData != null ? getNullableDate(sampleData.arrivalDateString()) : null;

        if (arrivalDate == null) {
            arrivalDate = preLimsArrivalDates.get(sample);
        }

        if (arrivalDate == null) {
            LOGGER.warn("Could not find a valid arrival date for sample: " + sample + " in LIMS");
        }

        return arrivalDate;
    }

    @Nullable
    public LocalDate samplingDate(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            final String samplingDateString = sampleData.samplingDateString();
            final LocalDate samplingDate = getNullableDate(samplingDateString);
            if (samplingDate == null && !samplesWithoutSamplingDate.contains(sample)) {
                LOGGER.warn("LIMS sampling date for " + sample + ": " + sampleData.samplingDateString() + " is not a valid date.");
            }
            return samplingDate;
        }
        return null;
    }

    @Nullable
    public Integer dnaNanograms(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
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
    public String purityShallowSeq(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        LimsShallowSeqData shallowSeq = shallowSeqPerSample.get(sample);

        if (sampleData != null) {
            boolean purityShallowExecuted = shallowSeqExecuted(sample);

            if (purityShallowExecuted && shallowSeq == null) {
                LOGGER.warn("BFX lims and lab status do not match for sample " + sample + "!");
            } else {
                if (purityShallowExecuted) {
                    LOGGER.info("Retrieved purity from shallow seq: " + shallowSeq.purityShallowSeq());
                    if (shallowSeq.purityShallowSeq().equals("below detection threshold")) {
                        return "below detection threshold";
                    } else {
                        try {
                            return Math.round(Double.parseDouble(shallowSeq.purityShallowSeq()) * 100) + "%";
                        } catch (final NumberFormatException e) {
                            LOGGER.warn("Could not convert shallow seq to a percentage: " + shallowSeq.purityShallowSeq());
                            return "N/A";
                        }
                    }
                } else {
                    return "not determined";
                }
            }
        }
        return "N/A";
    }

    @NotNull
    public String pathologyTumorPercentage(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            if (shallowSeqExecuted(sample)) {
                return "not determined";
            }

            String tumorPercentageString = sampleData.tumorPercentageString();
            if (tumorPercentageString == null) {
                return "N/A";
            } else {
                try {
                    return Long.toString(Math.round(Double.parseDouble(tumorPercentageString))) + "%";
                } catch (final NumberFormatException e) {
                    return "N/A";
                }
            }
        }
        return "N/A";
    }

    @NotNull
    public String primaryTumor(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            return sampleData.primaryTumor();
        }
        // No warning raised since initially this information was not tracked so this will be missing for early samples.
        return "N/A";
    }

    @NotNull
    public String labProcedures(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            return sampleData.labProcedures();
        }
        LOGGER.warn("Could not find lab SOP versions for sample: " + sample + " in LIMS");
        return "N/A";
    }

    @Nullable
    public LimsInformedConsent germlineFindigsWIDE(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            String germlineOptions = sampleData.germlineFindings();
            return LimsInformedConsent.extractChoiceInformedConsent(germlineOptions, sample);
        } else {
            return null;
        }
    }

    @NotNull
    public String hospitalPaSampleIdWIDE(@NotNull String sample){
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            return sampleData.hospitalPaSampleId();
        }
        LOGGER.info("hospital PA sample ID is not known");
        return "not known";
    }

    @Nullable
    private String submission(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.submission() : null;
    }

    private boolean shallowSeqExecuted(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        assert sampleData != null;

        String labRemarks = sampleData.labRemarks();

        if (sampleData.shallowSeq() == 1) {
            return true;
        } else if (labRemarks != null) {
            String labRemarksLowerCase = labRemarks.toLowerCase();
            return labRemarksLowerCase.contains("cpctwide") || labRemarksLowerCase.contains("shallowseq");
        }

        return false;
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
