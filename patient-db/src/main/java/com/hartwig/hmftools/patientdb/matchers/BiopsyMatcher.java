package com.hartwig.hmftools.patientdb.matchers;

import java.time.Duration;
import java.time.LocalDate;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.Config;
import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BiopsyMatcher {
    private static final Logger LOGGER = LogManager.getLogger(BiopsyMatcher.class);

    @NotNull
    public static List<BiopsyClinicalData> matchBiopsies(@NotNull final String patientId,
            @NotNull final List<BiopsyLimsData> sequencedBiopsies,
            @NotNull final List<BiopsyClinicalData> clinicalBiopsies) {
        if (clinicalBiopsies.size() < sequencedBiopsies.size()) {
            LOGGER.warn(patientId + ": contains less biopsies in ecrf (" + clinicalBiopsies.size()
                    + ") than biopsies sequenced (" + sequencedBiopsies.size() + ").");
        }
        final List<BiopsyClinicalData> matchedBiopsies = Lists.newArrayList();
        int startIndex = 0;
        for (final BiopsyLimsData sequencedBiopsy : sequencedBiopsies) {
            final MatchIndices indices = findMatchIndices(sequencedBiopsy, clinicalBiopsies, startIndex);
            if (indices.matchCount == 0) {
                LOGGER.warn(patientId + ": Could not match any clinical biopsy with sequenced biopsy: "
                        + sequencedBiopsy.sampleId() + "(" + sequencedBiopsy.samplingDate() + ","
                        + sequencedBiopsy.arrivalDate() + "): " + clinicalBiopsies.stream().map(
                        BiopsyClinicalData::date).collect(Collectors.toList()) + " on " + getMatchDateCriteria(
                        sequencedBiopsy));
                return clinicalBiopsies;
            } else if (indices.matchCount > 1) {
                LOGGER.warn(patientId + ": Found more than 1 possible clinical biopsy match for sequenced biopsy: "
                        + sequencedBiopsy.sampleId() + "(" + sequencedBiopsy.samplingDate() + ","
                        + sequencedBiopsy.arrivalDate() + "): " + clinicalBiopsies.stream().map(
                        BiopsyClinicalData::date).collect(Collectors.toList()) + " on " + getMatchDateCriteria(
                        sequencedBiopsy));
                return clinicalBiopsies;
            } else {
                for (int index = startIndex; index < indices.nextStartIndex; index++) {
                    final BiopsyClinicalData clinicalBiopsy = clinicalBiopsies.get(index);
                    if (index == indices.lastMatchIndex) {
                        matchedBiopsies.add(new BiopsyClinicalData(clinicalBiopsy.id(), clinicalBiopsy.date(),
                                clinicalBiopsy.location(), sequencedBiopsy.sampleId()));
                    } else
                        matchedBiopsies.add(clinicalBiopsy);
                }
            }
            startIndex = indices.nextStartIndex;
        }
        for (int index = startIndex; index < clinicalBiopsies.size(); index++) {
            matchedBiopsies.add(clinicalBiopsies.get(index));
        }
        return matchedBiopsies;
    }

    private static MatchIndices findMatchIndices(@NotNull final BiopsyLimsData sequencedBiopsy,
            @NotNull final List<BiopsyClinicalData> clinicalBiopsies, int startIndex) {
        int count = 0;
        int matchIndex = -1;
        int nextStartIndex = clinicalBiopsies.size();
        for (int index = startIndex; index < clinicalBiopsies.size(); index++) {
            final BiopsyClinicalData clinicalBiopsy = clinicalBiopsies.get(index);
            final LocalDate clinicalBiopsyDate = clinicalBiopsy.date();
            if (clinicalBiopsyDate == null) {
                count++;
            } else if (isWithinThreshold(sequencedBiopsy, clinicalBiopsy)) {
                count++;
                matchIndex = index;
            } else if (clinicalBiopsyDate.isAfter(sequencedBiopsy.date())) {
                nextStartIndex = index;
                break;
            }
        }
        return new MatchIndices(count, matchIndex, nextStartIndex);
    }

    private static class MatchIndices {
        final int matchCount;
        final int lastMatchIndex;
        final int nextStartIndex;

        MatchIndices(final int matchCount, final int lastMatchIndex, final int nextStartIndex) {
            this.matchCount = matchCount;
            this.lastMatchIndex = lastMatchIndex;
            this.nextStartIndex = nextStartIndex;
        }
    }

    private static boolean isWithinThreshold(@NotNull final BiopsyLimsData sequencedBiopsy,
            @NotNull final BiopsyClinicalData clinicalBiopsy) {
        final LocalDate biopsyDate = clinicalBiopsy.date();
        if (biopsyDate != null && (biopsyDate.isBefore(sequencedBiopsy.date()) || biopsyDate.isEqual(
                sequencedBiopsy.date()))) {
            final LocalDate limsSamplingDate = sequencedBiopsy.samplingDate();
            if (limsSamplingDate != null) {
                return Duration.between(biopsyDate.atStartOfDay(), limsSamplingDate.atStartOfDay()).toDays()
                        < Config.samplingDateThreshold;
            } else {
                return Duration.between(biopsyDate.atStartOfDay(),
                        sequencedBiopsy.arrivalDate().atStartOfDay()).toDays() < Config.arrivalDateThreshold;
            }
        }
        return false;
    }

    @NotNull
    private static String getMatchDateCriteria(@NotNull final BiopsyLimsData sequencedBiopsy) {
        if (sequencedBiopsy.samplingDate() != null) {
            return "sampling date " + sequencedBiopsy.samplingDate() + " threshold: " + Config.samplingDateThreshold;
        }
        return "arrival date " + sequencedBiopsy.arrivalDate() + " threshold: " + Config.arrivalDateThreshold;
    }
}
