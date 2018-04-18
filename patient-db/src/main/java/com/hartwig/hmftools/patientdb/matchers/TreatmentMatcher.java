package com.hartwig.hmftools.patientdb.matchers;

import java.time.Duration;
import java.time.LocalDate;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.Config;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class TreatmentMatcher {

    private TreatmentMatcher() {
    }

    @NotNull
    public static MatchResult<BiopsyTreatmentData> matchTreatmentsToBiopsies(@NotNull final String patientIdentifier,
            @NotNull final List<BiopsyData> biopsies, @NotNull final List<BiopsyTreatmentData> treatments) {
        final List<BiopsyTreatmentData> matchedTreatments = Lists.newArrayList();
        final List<ValidationFinding> findings = Lists.newArrayList();

        List<BiopsyData> remainingBiopsies = Lists.newArrayList(biopsies);

        Collections.sort(biopsies);
        Collections.sort(treatments);

        if (firstTreatmentBeforeFirstBiopsy(biopsies, treatments)) {
            findings.add(treatmentMatchFinding(patientIdentifier,
                    "First treatment prior to first biopsy",
                    "biopsy date: " + biopsies.get(0).date() + ", first treatment start: " + treatments.get(0).startDate()));
        }

        List<BiopsyTreatmentData> yesTreatments = getYesTreatments(treatments);
        List<BiopsyTreatmentData> notYesTreatments = getNotYesTreatments(treatments);

        // KODU: First match yes-treatments
        for (final BiopsyTreatmentData treatment : yesTreatments) {
            LocalDate startDate = treatment.startDate();
            if (startDate == null) {
                matchedTreatments.add(treatment);
            } else {
                BiopsyData bestMatch = null;
                for (BiopsyData remainingBiopsy : remainingBiopsies) {
                    if (isPossibleMatch(remainingBiopsy, startDate)) {
                        bestMatch = determineBestMatch(remainingBiopsy, bestMatch);
                    }
                }

                if (bestMatch != null) {
                    matchedTreatments.add(ImmutableBiopsyTreatmentData.builder().from(treatment).biopsyId(bestMatch.id()).build());
                    remainingBiopsies.remove(bestMatch);
                } else {
                    findings.add(treatmentMatchFinding(patientIdentifier,
                            "Could not find a biopsy match for a given treatment and having a start date!",
                            treatment.toString()));
                    matchedTreatments.add(treatment);
                }
            }
        }

        // KODU: Then distribute not-yes treatments over remaining biopsies.
        for (final BiopsyTreatmentData treatment : notYesTreatments) {
            final String treatmentGiven = treatment.treatmentGiven();

            if (treatmentGiven == null) {
                matchedTreatments.add(treatment);
            } else if (treatmentGiven.equalsIgnoreCase("no")) {
                if (remainingBiopsies.size() > 0) {
                    matchedTreatments.add(ImmutableBiopsyTreatmentData.builder()
                            .from(treatment)
                            .biopsyId(remainingBiopsies.get(0).id())
                            .build());
                    remainingBiopsies.remove(remainingBiopsies.get(0));
                } else {
                    matchedTreatments.add(treatment);
                }
            } else {
                matchedTreatments.add(treatment);
            }
        }
        

        return new MatchResult<>(matchedTreatments, findings);
    }

    private static boolean firstTreatmentBeforeFirstBiopsy(@NotNull final List<BiopsyData> sortedBiopsies,
            @NotNull final List<BiopsyTreatmentData> sortedTreatments) {
        if (sortedBiopsies.size() > 0 && sortedTreatments.size() > 0) {
            LocalDate firstBiopsyDate = sortedBiopsies.get(0).date();
            LocalDate firstTreatmentStart = sortedTreatments.get(0).startDate();

            if (firstBiopsyDate != null && firstTreatmentStart != null) {
                return firstTreatmentStart.isBefore(firstBiopsyDate);
            }
        }

        return false;
    }

    private static boolean isPossibleMatch(@NotNull final BiopsyData biopsy, @NotNull final LocalDate treatmentStartDate) {
        LocalDate biopsyDate = biopsy.date();

        return biopsyDate != null && (treatmentStartDate.isAfter(biopsyDate) || treatmentStartDate.isEqual(biopsyDate))
                && Duration.between(biopsyDate.atStartOfDay(), treatmentStartDate.atStartOfDay()).toDays()
                < Config.MAX_DAYS_BETWEEN_TREATMENT_AND_BIOPSY;
    }

    @NotNull
    private static BiopsyData determineBestMatch(@NotNull BiopsyData potentialBest, @Nullable BiopsyData currentBest) {
        if (currentBest == null) {
            return potentialBest;
        }

        LocalDate potentialBestBiopsyDate = potentialBest.date();
        LocalDate currentBestBiopsyDate = currentBest.date();

        assert potentialBestBiopsyDate != null && currentBestBiopsyDate != null;

        return potentialBestBiopsyDate.isAfter(currentBestBiopsyDate) ? potentialBest : currentBest;
    }

    @NotNull
    private static List<BiopsyTreatmentData> getYesTreatments(@NotNull List<BiopsyTreatmentData> treatments) {
        List<BiopsyTreatmentData> yesTreatments = Lists.newArrayList();
        for (BiopsyTreatmentData treatment : treatments) {
            String treatmentGiven = treatment.treatmentGiven();
            if (treatmentGiven != null && treatmentGiven.equalsIgnoreCase("yes")) {
                yesTreatments.add(treatment);
            }
        }
        return yesTreatments;
    }

    @NotNull
    private static List<BiopsyTreatmentData> getNotYesTreatments(@NotNull List<BiopsyTreatmentData> treatments) {
        List<BiopsyTreatmentData> notYesTreatments = Lists.newArrayList();
        for (BiopsyTreatmentData treatment : treatments) {
            String treatmentGiven = treatment.treatmentGiven();
            if (treatmentGiven == null || !treatmentGiven.equalsIgnoreCase("yes")) {
                notYesTreatments.add(treatment);
            }
        }
        return notYesTreatments;

    }

    @NotNull
    private static ValidationFinding treatmentMatchFinding(@NotNull String patientIdentifier, @NotNull String message,
            @NotNull String details) {
        return ValidationFinding.of("match", patientIdentifier, message, FormStatus.undefined(), details);
    }
}
