package com.hartwig.hmftools.patientdb.matchers;

import static com.hartwig.hmftools.patientdb.readers.BiopsyReader.FORM_BIOPS;

import java.time.Duration;
import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ImmutableValidationFinding;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
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
        if (clinicalBiopsies.size() < sequencedBiopsies.size()) {
            findings.add(new ImmutableValidationFinding("match", patientId, FORM_BIOPS,
                    "less biopsies in ecrf (" + clinicalBiopsies.size() + ") than biopsies sequenced (" + sequencedBiopsies.size() + ").",
                    "", ""));
        }
        List<BiopsyData> remainingBiopsies = clinicalBiopsies;
        for (final SampleData sequencedBiopsy : sequencedBiopsies) {
            final Map<Boolean, List<BiopsyData>> partitions = remainingBiopsies.stream()
                    .collect(Collectors.partitioningBy(clinicalBiopsy -> isPossibleMatch(sequencedBiopsy, clinicalBiopsy)));
            final List<BiopsyData> possibleMatches = partitions.get(true);
            if (possibleMatches.size() == 1 && possibleMatches.get(0).date() != null) {
                final BiopsyData clinicalBiopsy = partitions.get(true).get(0);
                matchedBiopsies.add(ImmutableBiopsyData.of(clinicalBiopsy.id(), clinicalBiopsy.date(), clinicalBiopsy.location(),
                        sequencedBiopsy.sampleId(), clinicalBiopsy.formStatus(), clinicalBiopsy.formLocked()));
                remainingBiopsies = partitions.get(false);
            } else if (possibleMatches.size() == 0 || (possibleMatches.size() == 1 && possibleMatches.get(0).date() == null)) {
                findings.add(new ImmutableValidationFinding("match", patientId, FORM_BIOPS,
                        "could not match any clinical biopsy with sequenced biopsy: " + sequencedBiopsy.sampleId() + "("
                                + sequencedBiopsy.samplingDate() + "," + sequencedBiopsy.arrivalDate() + "): " + clinicalBiopsies.stream()
                                .map(BiopsyData::date)
                                .collect(Collectors.toList()) + " on " + getMatchDateCriteria(sequencedBiopsy), "", ""));
                // MIVO: abort finding new matches if we can't match one sequenced biopsy
                return new MatchResult<>(clinicalBiopsies, findings);
            } else if (possibleMatches.size() > 1) {
                findings.add(new ImmutableValidationFinding("match", patientId, FORM_BIOPS,
                        "more than 1 possible clinical biopsy match for sequenced biopsy: " + sequencedBiopsy.sampleId() + "("
                                + sequencedBiopsy.samplingDate() + "," + sequencedBiopsy.arrivalDate() + "): " + clinicalBiopsies.stream()
                                .map(BiopsyData::date)
                                .collect(Collectors.toList()) + " on " + getMatchDateCriteria(sequencedBiopsy), "", ""));
                // MIVO: abort finding new matches if we can't match one sequenced biopsy
                return new MatchResult<>(clinicalBiopsies, findings);
            }
        }
        matchedBiopsies.addAll(remainingBiopsies);
        return new MatchResult<>(matchedBiopsies, findings);
    }

    private static boolean isPossibleMatch(@NotNull final SampleData sequencedBiopsy, @NotNull final BiopsyData clinicalBiopsy) {
        return clinicalBiopsy.date() == null || isWithinThreshold(sequencedBiopsy, clinicalBiopsy);
    }

    private static boolean isWithinThreshold(@NotNull final SampleData sequencedBiopsy, @NotNull final BiopsyData clinicalBiopsy) {
        final LocalDate biopsyDate = clinicalBiopsy.date();
        if (biopsyDate != null && (biopsyDate.isBefore(sequencedBiopsy.date()) || biopsyDate.isEqual(sequencedBiopsy.date()))) {
            final LocalDate limsSamplingDate = sequencedBiopsy.samplingDate();
            if (limsSamplingDate != null) {
                return Duration.between(biopsyDate.atStartOfDay(), limsSamplingDate.atStartOfDay()).toDays()
                        < Config.SAMPLING_DATE_THRESHOLD;
            } else {
                return Duration.between(biopsyDate.atStartOfDay(), sequencedBiopsy.arrivalDate().atStartOfDay()).toDays()
                        < Config.ARRIVAL_DATE_THRESHOLD;
            }
        }
        return false;
    }

    @NotNull
    private static String getMatchDateCriteria(@NotNull final SampleData sequencedBiopsy) {
        if (sequencedBiopsy.samplingDate() != null) {
            return "sampling date " + sequencedBiopsy.samplingDate() + " threshold: " + Config.SAMPLING_DATE_THRESHOLD;
        }
        return "arrival date " + sequencedBiopsy.arrivalDate() + " threshold: " + Config.ARRIVAL_DATE_THRESHOLD;
    }
}
