package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.amber.AmberSample.DO_NOT_MATCH;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class AmberPatientFactory {

    public static List<AmberPatient> create(double matchCutoff, @NotNull final AmberSample victim,
            @NotNull final List<AmberSample> existingSamples, @NotNull final List<AmberPatient> existingPatients) {
        final List<AmberPatient> result = Lists.newArrayList();
        final String sampleId = victim.sampleId();

        final List<SiteMatch> matches = allMatches(victim, existingSamples);
        final Set<String> patientSamples = Sets.newHashSet();
        for (SiteMatch match : matches) {
            if (Doubles.greaterOrEqual(match.likelihood(), matchCutoff)) {
                patientSamples.add(match.victim);
                patientSamples.add(match.other);
            }
        }

        if (!patientSamples.isEmpty()) {
            final List<AmberPatient> matchedPatients = existingPatients.stream()
                    .filter(x -> patientSamples.contains(x.sample()) || patientSamples.contains(x.otherSample()))
                    .collect(Collectors.toList());

            matchedPatients.forEach(x -> {
                patientSamples.add(x.sample());
                patientSamples.add(x.otherSample());
            });

            final Set<Integer> matchedPatientId = matchedPatients.stream().map(AmberPatient::patientId).collect(Collectors.toSet());
            final int patientId = matchedPatientId.size() == 1
                    ? matchedPatientId.iterator().next()
                    : existingPatients.stream().mapToInt(AmberPatient::patientId).max().orElse(0) + 1;

            // All existing matches (excluding victim)
            matchedPatients.stream()
                    .filter(x -> !x.sample().equals(sampleId) && !x.otherSample().equals(sampleId))
                    .map(x -> ImmutableAmberPatient.builder().from(x).patientId(patientId).build())
                    .forEach(result::add);

            // Add new matches
            matches.stream()
                    .filter(x -> patientSamples.contains(x.other))
                    .map(x -> create(patientId, x))
                    .forEach(result::add);

        }

        return result;
    }

    @NotNull
    private static AmberPatient create(int patientId, @NotNull final SiteMatch match) {
        final String first = match.victim.compareTo(match.other) < 0 ? match.victim : match.other;
        final String second = first.equals(match.victim) ? match.other : match.victim;

        return ImmutableAmberPatient.builder()
                .patientId(patientId)
                .sample(first)
                .otherSample(second)
                .sites(match.sites)
                .matches(match.matches)
                .build();
    }

    @NotNull
    private static List<SiteMatch> allMatches(@NotNull final AmberSample victim, @NotNull final List<AmberSample> existing) {
        List<SiteMatch> result = Lists.newArrayList();
        for (AmberSample existingSample : existing) {
            if (!victim.sampleId().equals(existingSample.sampleId())) {
                result.add(create(victim, existingSample));
            }
        }
        return new ArrayList<>(result);
    }

    @NotNull
    private static SiteMatch create(@NotNull final AmberSample victim, @NotNull final AmberSample other) {
        final byte[] entries = victim.entries();
        byte[] otherEntries = other.entries();
        if (entries.length != otherEntries.length) {
            throw new IllegalArgumentException("Unable to match different sized identities");
        }

        int matches = 0;
        int sites = 0;
        for (int i = 0; i < entries.length; i++) {
            byte myByte = entries[i];
            byte otherByte = otherEntries[i];

            if (myByte != DO_NOT_MATCH && otherByte != DO_NOT_MATCH) {
                sites++;
                matches += (myByte == otherByte ? 1 : 0);
            }
        }

        return new SiteMatch(victim.sampleId(), other.sampleId(), matches, sites);
    }

    private static class SiteMatch {
        private final String victim;
        private final String other;
        private final int matches;
        private final int sites;

        SiteMatch(final String firstSampleId, final String secondSampleId, final int matches, final int sites) {
            this.victim = firstSampleId;
            this.other = secondSampleId;
            this.matches = matches;
            this.sites = sites;
        }

        public double likelihood() {
            return matches / (double) sites;
        }
    }

}
