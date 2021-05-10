package com.hartwig.hmftools.patientdb.clinical.matchers;

import java.time.Duration;
import java.time.LocalDate;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.clinical.ClinicalConstants;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableBiopsyTreatmentData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class TreatmentMatcher {

    private TreatmentMatcher() {
    }

    @NotNull
    public static MatchResult<BiopsyTreatmentData> matchTreatmentsToBiopsies(@NotNull String patientIdentifier,
            @NotNull List<BiopsyData> biopsies, @NotNull List<BiopsyTreatmentData> treatments) {
        List<BiopsyTreatmentData> matchedTreatments = Lists.newArrayList();
        List<ValidationFinding> findings = Lists.newArrayList();

        List<BiopsyData> remainingBiopsies = Lists.newArrayList(biopsies);

        Collections.sort(biopsies);
        Collections.sort(treatments);

        List<BiopsyTreatmentData> yesTreatments = getYesTreatments(treatments);
        List<BiopsyTreatmentData> notYesTreatments = getNotYesTreatments(treatments);

        // First match yes-treatments
        for (BiopsyTreatmentData treatment : yesTreatments) {
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
                    matchedTreatments.add(treatment);
                }
            }
        }

        // Then distribute not-yes treatments over remaining biopsies.
        for (BiopsyTreatmentData treatment : notYesTreatments) {
            String treatmentGiven = treatment.treatmentGiven();

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

        findings.addAll(validateMatchingForMatchedBiopsies(patientIdentifier, matchedTreatments, biopsies));

        return new MatchResult<>(matchedTreatments, findings);
    }

    @NotNull
    private static Collection<ValidationFinding> validateMatchingForMatchedBiopsies(@NotNull String patientIdentifier,
            @NotNull List<BiopsyTreatmentData> matchedTreatments, @NotNull List<BiopsyData> biopsies) {
        List<ValidationFinding> findings = Lists.newArrayList();
        List<Integer> matchedBiopsyIds = Lists.newArrayList();

        for (BiopsyTreatmentData treatment : matchedTreatments) {
            if (treatment.biopsyId() != null) {
                matchedBiopsyIds.add(treatment.biopsyId());
            }
        }

        for (BiopsyData biopsy : biopsies) {
            if (biopsy.sampleId() != null && !matchedBiopsyIds.contains(biopsy.id())) {
                findings.add(treatmentMatchFinding(patientIdentifier, "Could not match treatment to matched biopsy", biopsy.toString()));
            }
        }
        return findings;
    }

    private static boolean isPossibleMatch(@NotNull BiopsyData biopsy, @NotNull LocalDate treatmentStartDate) {
        LocalDate biopsyDate = biopsy.date();

        return biopsyDate != null && (treatmentStartDate.isAfter(biopsyDate) || treatmentStartDate.isEqual(biopsyDate))
                && Duration.between(biopsyDate.atStartOfDay(), treatmentStartDate.atStartOfDay()).toDays()
                < ClinicalConstants.MAX_DAYS_BETWEEN_TREATMENT_AND_BIOPSY;
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
