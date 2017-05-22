package com.hartwig.hmftools.patientdb.matchers;

import java.time.Duration;
import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.Config;
import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class TreatmentMatcher {
    private static final Logger LOGGER = LogManager.getLogger(TreatmentMatcher.class);

    private TreatmentMatcher() {
    }

    @NotNull
    public static List<BiopsyTreatmentData> matchTreatments(@NotNull final String patientId,
            @NotNull final List<BiopsyClinicalData> biopsies, @NotNull final List<BiopsyTreatmentData> treatments) {
        final List<BiopsyTreatmentData> matchedTreatments = Lists.newArrayList();
        if (biopsies.size() < treatments.size()) {
            LOGGER.warn(patientId + ": has more treatments(" + treatments.size() + ") than biopsies(" + biopsies.size()
                    + ").");
        }
        List<BiopsyClinicalData> remainingBiopsies = biopsies;
        for (final BiopsyTreatmentData treatment : treatments) {
            final String treatmentGiven = treatment.treatmentGiven();
            if (treatmentGiven == null || treatmentGiven.toLowerCase().equals("no")) {
                matchedTreatments.add(treatment);
                continue;
            }
            final Map<Boolean, List<BiopsyClinicalData>> partitions = remainingBiopsies.stream().collect(
                    Collectors.partitioningBy(
                            clinicalBiopsy -> isPossibleMatch(clinicalBiopsy.date(), treatment.startDate())));
            final List<BiopsyClinicalData> possibleMatches = partitions.get(true);
            if (possibleMatches.size() == 0) {
                LOGGER.warn(patientId + ": no biopsy match for treatment [" + treatment.startDate() + " - "
                        + treatment.endDate() + "]");
                matchedTreatments.add(treatment);
            } else if (possibleMatches.size() > 1) {
                LOGGER.warn(patientId + ": multiple biopsy matches for treatment [" + treatment.startDate() + " - "
                        + treatment.endDate() + "]: " + possibleMatches.stream().map(
                        biopsy -> "" + biopsy.date()).collect(Collectors.toList()));
                matchedTreatments.add(treatment);
            } else if (possibleMatches.size() == 1 && possibleMatches.get(0).date() == null) {
                LOGGER.warn(patientId + ": treatment [" + treatment.startDate() + " - " + treatment.endDate()
                        + "] matched biopsy with null date.");
                matchedTreatments.add(treatment);
            } else {
                final BiopsyClinicalData clinicalBiopsy = possibleMatches.get(0);
                matchedTreatments.add(
                        new BiopsyTreatmentData(treatment.id(), treatment.treatmentGiven(), treatment.startDate(),
                                treatment.endDate(), treatment.drugs(), clinicalBiopsy.id()));
                remainingBiopsies = partitions.get(false);
            }
        }
        return matchedTreatments;
    }

    private static boolean isPossibleMatch(@Nullable final LocalDate biopsyDate,
            @Nullable final LocalDate treatmentStartDate) {
        return biopsyDate == null || isWithinThreshold(biopsyDate, treatmentStartDate);
    }

    private static boolean isWithinThreshold(@Nullable final LocalDate biopsyDate,
            @Nullable final LocalDate treatmentStartDate) {
        return biopsyDate != null && treatmentStartDate != null && (treatmentStartDate.isAfter(biopsyDate)
                || treatmentStartDate.isEqual(biopsyDate))
                && Duration.between(biopsyDate.atStartOfDay(), treatmentStartDate.atStartOfDay()).toDays()
                < Config.maxDaysBetweenTreatmentAndBiopsy;
    }
}
