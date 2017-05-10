package com.hartwig.hmftools.patientdb.matchers;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TreatmentResponseMatcher {

    private static final Logger LOGGER = LogManager.getLogger(TreatmentResponseMatcher.class);

    public static List<BiopsyTreatmentResponseData> matchTreatmentResponses(@NotNull final String patientId,
            @NotNull final List<BiopsyTreatmentData> treatments,
            @NotNull final List<BiopsyTreatmentResponseData> responses) {
        final List<BiopsyTreatmentResponseData> matchedResponses = Lists.newArrayList();
        if (treatments.size() == 1) {
            for (final BiopsyTreatmentResponseData response : responses) {
                matchedResponses.add(new BiopsyTreatmentResponseData(treatments.get(0).id(), response.assessmentDate(),
                        response.responseDate(), response.response(), response.measurementDone()));
            }
        } else {
            for (final BiopsyTreatmentResponseData response : responses) {
                final Integer treatmentId = findTreatmentIdForResponse(patientId, response, treatments);
                matchedResponses.add(new BiopsyTreatmentResponseData(treatmentId, response.assessmentDate(),
                        response.responseDate(), response.response(), response.measurementDone()));
            }
        }
        return matchedResponses;
    }

    @Nullable
    private static Integer findTreatmentIdForResponse(@NotNull final String patientId,
            @NotNull final BiopsyTreatmentResponseData response, @NotNull final List<BiopsyTreatmentData> treatments) {
        BiopsyTreatmentData matchedTreatment = null;
        for (final BiopsyTreatmentData treatment : treatments) {
            if (responseMatchesTreatment(response, treatment)) {
                if (matchedTreatment == null) {
                    matchedTreatment = treatment;
                } else {
                    LOGGER.warn(patientId + ": treatment response(" + response.date() + ") matched with: "
                            + matchedTreatment.id() + " (" + matchedTreatment.startDate() + ","
                            + matchedTreatment.endDate() + ") can also be matched with: " + treatment.id() + " ("
                            + treatment.startDate() + "," + treatment.endDate() + ")");
                }
            }
        }
        if (matchedTreatment != null) {
            return matchedTreatment.id();
        }
        return null;
    }

    private static boolean responseMatchesTreatment(@NotNull final BiopsyTreatmentResponseData response,
            @NotNull final BiopsyTreatmentData treatment) {
        final LocalDate treatmentStart = treatment.startDate();
        final LocalDate treatmentEnd = treatment.endDate();
        final LocalDate responseDate = response.date();
        if (treatmentStart == null || responseDate == null) {
            return false;
        }
        return (responseDate.isAfter(treatmentStart) && (treatmentEnd == null || responseDate.isBefore(treatmentEnd)
                || responseDate.isEqual(treatmentEnd)));
    }
}
