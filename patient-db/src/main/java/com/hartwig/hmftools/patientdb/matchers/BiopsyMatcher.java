package com.hartwig.hmftools.patientdb.matchers;

import java.time.Duration;
import java.time.LocalDate;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BiopsyMatcher {
    private static final Logger LOGGER = LogManager.getLogger(BiopsyMatcher.class);

    @NotNull
    public static List<BiopsyClinicalData> matchBiopsies(@NotNull final String patientId,
            @NotNull final List<BiopsyLimsData> sequencedBiopsies,
            @NotNull final List<BiopsyClinicalData> clinicalBiopsies) {
        if (clinicalBiopsies.size() == sequencedBiopsies.size()) {
            // MIVO: if there's n biopsies sequenced, and n (non-empty) entries in ecrf => match
            return matchBiopsiesByIndex(sequencedBiopsies, clinicalBiopsies);
        } else {
            if (clinicalBiopsies.size() < sequencedBiopsies.size()) {
                LOGGER.warn("Patient: " + patientId + " contains less entries in ecrf (" + clinicalBiopsies.size()
                        + ") than biopsies sequenced (" + sequencedBiopsies.size() + ").");
            }
            return matchBiopsiesByDate(sequencedBiopsies, clinicalBiopsies);
        }
    }

    @NotNull
    private static List<BiopsyClinicalData> matchBiopsiesByIndex(@NotNull final List<BiopsyLimsData> sequencedBiopsies,
            @NotNull final List<BiopsyClinicalData> clinicalBiopsies) {
        final List<BiopsyClinicalData> matchedBiopsies = Lists.newArrayList();
        for (int index = 0; index < clinicalBiopsies.size(); index++) {
            final BiopsyClinicalData clinicalBiopsy = clinicalBiopsies.get(index);
            final String sampleId = sequencedBiopsies.get(index).sampleId();
            matchedBiopsies.add(
                    new BiopsyClinicalData(clinicalBiopsy.id(), clinicalBiopsy.date(), clinicalBiopsy.location(),
                            clinicalBiopsy.treatment(), sampleId));
        }
        return matchedBiopsies;
    }

    @NotNull
    private static List<BiopsyClinicalData> matchBiopsiesByDate(@NotNull final List<BiopsyLimsData> sequencedBiopsies,
            @NotNull final List<BiopsyClinicalData> clinicalBiopsies) {
        final List<BiopsyClinicalData> matchedBiopsies = Lists.newArrayList();
        for (final BiopsyLimsData biopsyLimsData : sequencedBiopsies) {
            final LocalDate biopsyArrivalDate = biopsyLimsData.arrivalDate();
            BiopsyClinicalData sequencedBiopsyClinicalData = findBiopsyBeforeArrivalDate(biopsyArrivalDate,
                    clinicalBiopsies);
            if (sequencedBiopsyClinicalData != null) {
                checkDurationBetweenDates(biopsyArrivalDate, sequencedBiopsyClinicalData.date());
                final List<Integer> matchedBiopsiesIds = matchedBiopsies.stream().map(BiopsyClinicalData::id).collect(
                        Collectors.toList());
                if (!matchedBiopsiesIds.contains(sequencedBiopsyClinicalData.id())) {
                    matchedBiopsies.add(new BiopsyClinicalData(sequencedBiopsyClinicalData.id(),
                            sequencedBiopsyClinicalData.date(), sequencedBiopsyClinicalData.location(),
                            sequencedBiopsyClinicalData.treatment(), biopsyLimsData.sampleId()));
                }
            }
        }
        final List<Integer> matchedBiopsiesIds = matchedBiopsies.stream().map(BiopsyClinicalData::id).collect(
                Collectors.toList());
        for (final BiopsyClinicalData biopsyClinicalData : clinicalBiopsies) {
            if (!matchedBiopsiesIds.contains(biopsyClinicalData.id())) {
                matchedBiopsies.add(biopsyClinicalData);
            }
        }
        return matchedBiopsies;
    }

    @Nullable
    private static BiopsyClinicalData findBiopsyBeforeArrivalDate(@NotNull final LocalDate arrivalDate,
            @NotNull final List<BiopsyClinicalData> clinicalBiopsies) {
        BiopsyClinicalData result = null;
        for (final BiopsyClinicalData biopsyClinicalData : clinicalBiopsies) {
            final LocalDate biopsyDate = biopsyClinicalData.date();
            if (biopsyDate != null && biopsyDate.isBefore(arrivalDate)) {
                result = biopsyClinicalData;
            }
        }
        return result;
    }

    private static void checkDurationBetweenDates(@NotNull final LocalDate arrivalDate,
            @Nullable final LocalDate biopsyDate) {
        final int maxDaysBetweenBiopsyDateAndArrival = 180;
        if (biopsyDate != null && Duration.between(biopsyDate.atStartOfDay(), arrivalDate.atStartOfDay()).toDays()
                > maxDaysBetweenBiopsyDateAndArrival) {
            LOGGER.warn("Time between biopsy date (" + biopsyDate + ") and biopsy arrival date (" + arrivalDate
                    + " is greater than " + maxDaysBetweenBiopsyDateAndArrival + " days.");
        }
    }
}
