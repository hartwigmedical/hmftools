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
    private final Map<String, LimsShallowSeqData> dataShallowSeq;

    Lims(@NotNull final Map<String, LimsJsonSampleData> dataPerSample, @NotNull final Map<String, LimsJsonSubmissionData> dataPerSubmission,
            @NotNull final Map<String, LocalDate> preLimsArrivalDates, @NotNull final Set<String> samplesWithoutSamplingDate,
            @NotNull final Map<String, LimsShallowSeqData> dataShallowSeq) {
        this.dataPerSample = dataPerSample;
        this.dataPerSubmission = dataPerSubmission;
        this.preLimsArrivalDates = preLimsArrivalDates;
        this.samplesWithoutSamplingDate = samplesWithoutSamplingDate;
        this.dataShallowSeq = dataShallowSeq;
    }

    public int sampleCount() {
        return dataPerSample.size();
    }

    @NotNull
    public String contactEmail(@NotNull final String sample) {
        String submission = submission(sample);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.contactEmail() : "N/A";
    }

    @NotNull
    public String contactName(@NotNull final String sample) {
        String submission = submission(sample);
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        return submissionData != null ? submissionData.contactName() : "N/A";
    }

    @Nullable
    public String patientNumber(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.patientNumber() : null;
    }

    @NotNull
    public String label(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.labelSample() : "N/A";
    }

    public boolean isCoreSample(@NotNull final String sample) {
        String label = label(sample);
        return label.equalsIgnoreCase("core");
    }

    @NotNull
    public String projectName(@NotNull final String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.projectName() : "N/A";
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
    public String purityShallowSeq(@NotNull String sample, boolean isSequenced) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        LimsShallowSeqData shallowSeq = dataShallowSeq.get(sample);

        if (sampleData != null) {
            boolean purityShallowExecuted = shallowSeqExecuted(sample);

            if (purityShallowExecuted && shallowSeq == null) {
                LOGGER.warn("BFX lims and lab status do not match for sample " + sample + "!");
            } else if (!isSequenced) {
                if (purityShallowExecuted) {
                    LOGGER.info("Retrieved purity from shallow seq.");
                    return Math.round(shallowSeq.purityShallowSeq() * 100) + "%";
                } else {
                    LOGGER.info("Pathology tumor percentage retrieved for sample as proxy for shallow seq purity.");
                    return pathologyTumorPercentage(sample);
                }
            } else {
                LOGGER.info("Pathology tumor percentage retrieved for sample as proxy for shallow seq purity.");
                return pathologyTumorPercentage(sample);

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
    private String submission(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        return sampleData != null ? sampleData.submission() : null;
    }

    private boolean shallowSeqExecuted(@NotNull String sample) {
        LimsJsonSampleData sampleData = dataPerSample.get(sample);
        assert sampleData != null;

        String labRemarks = sampleData.labRemarks();

        return isCoreSample(sample) || (labRemarks != null && (labRemarks.contains("CPCTWIDE") || labRemarks.contains("ShallowSeq")));
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
