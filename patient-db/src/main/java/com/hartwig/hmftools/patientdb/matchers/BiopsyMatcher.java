package com.hartwig.hmftools.patientdb.matchers;

import java.time.Duration;
import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.Config;
import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class BiopsyMatcher {
    private static final Logger LOGGER = LogManager.getLogger(BiopsyMatcher.class);

    private BiopsyMatcher() {
    }

    @NotNull
    public static List<BiopsyClinicalData> matchBiopsies(@NotNull final String patientId,
            @NotNull final List<BiopsyLimsData> sequencedBiopsies,
            @NotNull final List<BiopsyClinicalData> clinicalBiopsies) {
        final List<BiopsyClinicalData> matchedBiopsies = Lists.newArrayList();
        if (clinicalBiopsies.size() < sequencedBiopsies.size()) {
            LOGGER.warn(patientId + ": contains less biopsies in ecrf (" + clinicalBiopsies.size()
                    + ") than biopsies sequenced (" + sequencedBiopsies.size() + ").");
        }
        List<BiopsyClinicalData> remainingBiopsies = clinicalBiopsies;
        for (final BiopsyLimsData sequencedBiopsy : sequencedBiopsies) {
            final Map<Boolean, List<BiopsyClinicalData>> partitions = remainingBiopsies.stream().collect(
                    Collectors.partitioningBy(clinicalBiopsy -> isPossibleMatch(sequencedBiopsy, clinicalBiopsy)));
            final List<BiopsyClinicalData> possibleMatches = partitions.get(true);
            if (possibleMatches.size() == 1 && possibleMatches.get(0).date() != null) {
                final BiopsyClinicalData clinicalBiopsy = partitions.get(true).get(0);
                matchedBiopsies.add(
                        new BiopsyClinicalData(clinicalBiopsy.id(), clinicalBiopsy.date(), clinicalBiopsy.location(),
                                sequencedBiopsy.sampleId()));
                remainingBiopsies = partitions.get(false);
            } else if (possibleMatches.size() == 0 || (possibleMatches.size() == 1
                    && possibleMatches.get(0).date() == null)) {
                LOGGER.warn(patientId + ": Could not match any clinical biopsy with sequenced biopsy: "
                        + sequencedBiopsy.sampleId() + "(" + sequencedBiopsy.samplingDate() + ","
                        + sequencedBiopsy.arrivalDate() + "): " + clinicalBiopsies.stream().map(
                        BiopsyClinicalData::date).collect(Collectors.toList()) + " on " + getMatchDateCriteria(
                        sequencedBiopsy));
                // MIVO: abort finding new matches if we can't match one sequenced biopsy
                return clinicalBiopsies;
            } else if (possibleMatches.size() > 1) {
                LOGGER.warn(patientId + ": Found more than 1 possible clinical biopsy match for sequenced biopsy: "
                        + sequencedBiopsy.sampleId() + "(" + sequencedBiopsy.samplingDate() + ","
                        + sequencedBiopsy.arrivalDate() + "): " + clinicalBiopsies.stream().map(
                        BiopsyClinicalData::date).collect(Collectors.toList()) + " on " + getMatchDateCriteria(
                        sequencedBiopsy));
                // MIVO: abort finding new matches if we can't match one sequenced biopsy
                return clinicalBiopsies;
            }
        }
        matchedBiopsies.addAll(remainingBiopsies);
        return matchedBiopsies;
    }

    private static boolean isPossibleMatch(@NotNull final BiopsyLimsData sequencedBiopsy,
            @NotNull final BiopsyClinicalData clinicalBiopsy) {
        return clinicalBiopsy.date() == null || isWithinThreshold(sequencedBiopsy, clinicalBiopsy);
    }

    private static boolean isWithinThreshold(@NotNull final BiopsyLimsData sequencedBiopsy,
            @NotNull final BiopsyClinicalData clinicalBiopsy) {
        final LocalDate biopsyDate = clinicalBiopsy.date();
        if (biopsyDate != null && (biopsyDate.isBefore(sequencedBiopsy.date()) || biopsyDate.isEqual(
                sequencedBiopsy.date()))) {
            final LocalDate limsSamplingDate = sequencedBiopsy.samplingDate();
            if (limsSamplingDate != null) {
                return Duration.between(biopsyDate.atStartOfDay(), limsSamplingDate.atStartOfDay()).toDays()
                        < Config.SAMPLING_DATE_THRESHOLD;
            } else {
                return Duration.between(biopsyDate.atStartOfDay(),
                        sequencedBiopsy.arrivalDate().atStartOfDay()).toDays() < Config.ARRIVAL_DATE_THRESHOLD;
            }
        }
        return false;
    }

    @NotNull
    private static String getMatchDateCriteria(@NotNull final BiopsyLimsData sequencedBiopsy) {
        if (sequencedBiopsy.samplingDate() != null) {
            return "sampling date " + sequencedBiopsy.samplingDate() + " threshold: " + Config.SAMPLING_DATE_THRESHOLD;
        }
        return "arrival date " + sequencedBiopsy.arrivalDate() + " threshold: " + Config.ARRIVAL_DATE_THRESHOLD;
    }
}
