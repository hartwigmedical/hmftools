package com.hartwig.hmftools.patientdb.matchers;

import java.time.Duration;
import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.Config;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyData;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.jetbrains.annotations.NotNull;

public final class BiopsyMatcher {
    private BiopsyMatcher() {
    }

    @NotNull
    public static MatchResult<BiopsyData> matchBiopsiesToTumorSamples(@NotNull final String patientId,
            @NotNull final List<SampleData> sequencedBiopsies, @NotNull final List<BiopsyData> clinicalBiopsies) {
        final List<BiopsyData> matchedBiopsies = Lists.newArrayList();
        final List<ValidationFinding> findings = Lists.newArrayList();

        List<BiopsyData> remainingBiopsies = clinicalBiopsies;
        if (clinicalBiopsies.size() < sequencedBiopsies.size()) {
            findings.add(biopsyMatchFinding(patientId,
                    "Not enough clinical biopsy forms to match for every sequenced sample",
                    "Clinical biopsies: " + clinicalBiopsies.size() + ", sequenced samples: " + sequencedBiopsies.size()));
        }

        for (final SampleData sequencedBiopsy : sequencedBiopsies) {
            final Map<Boolean, List<BiopsyData>> partitions = remainingBiopsies.stream()
                    .collect(Collectors.partitioningBy(clinicalBiopsy -> isPossibleMatch(sequencedBiopsy, clinicalBiopsy)));
            final List<BiopsyData> possibleMatches = partitions.get(true);
            if (possibleMatches.size() == 1 && possibleMatches.get(0).date() != null) {
                final BiopsyData clinicalBiopsy = possibleMatches.get(0);
                matchedBiopsies.add(ImmutableBiopsyData.builder().from(clinicalBiopsy).sampleId(sequencedBiopsy.sampleId()).build());
                if (hasMismatchOnSamplingDate(sequencedBiopsy, clinicalBiopsy)) {
                    findings.add(biopsyMatchFinding(patientId,
                            "Sampling date does not equal biopsy date in matched biopsy",
                            "sampling date: " + sequencedBiopsy.samplingDate() + "; ecrf biopsy date: " + clinicalBiopsy.date()));
                }
                remainingBiopsies = partitions.get(false);
            } else if (possibleMatches.size() == 0 || (possibleMatches.size() == 1 && possibleMatches.get(0).date() == null)) {
                findings.add(biopsyMatchFinding(patientId,
                        "Could not match any clinical biopsy with sequenced sample.",
                        "sample: " + sampleDataToString(sequencedBiopsy) + "; ecrf biopsies: " + clinicalBiopsies.stream()
                                .map(BiopsyData::date)
                                .collect(Collectors.toList()) + ". match criteria: " + getMatchDateCriteria(sequencedBiopsy)));
                // MIVO: abort finding new matches if we can't match one sequenced biopsy
                return new MatchResult<>(clinicalBiopsies, findings);
            } else if (possibleMatches.size() > 1) {
                findings.add(biopsyMatchFinding(patientId,
                        "More than 1 possible clinical biopsy match for sequenced sample.",
                        "sample: " + sampleDataToString(sequencedBiopsy) + "; ecrf biopsies: " + clinicalBiopsies.stream()
                                .map(BiopsyData::date)
                                .collect(Collectors.toList()) + ". match criteria: " + getMatchDateCriteria(sequencedBiopsy)));
                // MIVO: abort finding new matches if we can't match one sequenced biopsy
                return new MatchResult<>(clinicalBiopsies, findings);
            }
        }
        matchedBiopsies.addAll(remainingBiopsies);
        return new MatchResult<>(matchedBiopsies, findings);
    }

    @NotNull
    private static String sampleDataToString(@NotNull final SampleData sampleData) {
        return sampleData.sampleId() + "(S:" + sampleData.samplingDate() + ", A:" + sampleData.arrivalDate() + ")";
    }

    private static boolean isPossibleMatch(@NotNull final SampleData sequencedBiopsy, @NotNull final BiopsyData clinicalBiopsy) {
        return clinicalBiopsy.date() == null || (isWithinThreshold(sequencedBiopsy, clinicalBiopsy)
                && clinicalBiopsy.isPotentiallyEvaluable());
    }

    private static boolean isWithinThreshold(@NotNull final SampleData sequencedBiopsy, @NotNull final BiopsyData clinicalBiopsy) {
        final LocalDate biopsyDate = clinicalBiopsy.date();
        assert biopsyDate != null;

        final LocalDate samplingDate = sequencedBiopsy.samplingDate();
        if (samplingDate != null) {
            return Math.abs(Duration.between(biopsyDate.atStartOfDay(), samplingDate.atStartOfDay()).toDays())
                    < Config.MAX_DAYS_BETWEEN_SAMPLING_AND_BIOPSY_DATE;
        } else {
            long daysFromBiopsyTakenToArrival =
                    Duration.between(biopsyDate.atStartOfDay(), sequencedBiopsy.arrivalDate().atStartOfDay()).toDays();
            return daysFromBiopsyTakenToArrival >= 0 && daysFromBiopsyTakenToArrival < Config.MAX_DAYS_ARRIVAL_DATE_AFTER_BIOPSY_DATE;
        }
    }

    private static boolean hasMismatchOnSamplingDate(@NotNull SampleData sequencedBiopsy, @NotNull BiopsyData clinicalBiopsy) {
        LocalDate samplingDate = sequencedBiopsy.samplingDate();
        LocalDate biopsyDate = clinicalBiopsy.date();

        if (samplingDate != null) {
            assert biopsyDate != null;
            return !biopsyDate.equals(samplingDate);
        }

        return false;
    }

    @NotNull
    private static String getMatchDateCriteria(@NotNull final SampleData sampleData) {
        if (sampleData.samplingDate() != null) {
            return "Sampling date " + sampleData.samplingDate() + " threshold: " + Config.MAX_DAYS_BETWEEN_SAMPLING_AND_BIOPSY_DATE;
        }
        return "Arrival date " + sampleData.arrivalDate() + " threshold: " + Config.MAX_DAYS_ARRIVAL_DATE_AFTER_BIOPSY_DATE;
    }

    @NotNull
    private static ValidationFinding biopsyMatchFinding(@NotNull String patientIdentifier, @NotNull String finding,
            @NotNull String details) {
        return ValidationFinding.of("match", patientIdentifier, finding, FormStatus.undefined(), details);
    }
}
