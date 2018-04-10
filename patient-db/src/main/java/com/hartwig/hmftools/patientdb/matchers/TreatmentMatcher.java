package com.hartwig.hmftools.patientdb.matchers;

import static com.hartwig.hmftools.patientdb.readers.cpct.BiopsyTreatmentReader.FORM_TREATMENT;

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
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class TreatmentMatcher {

    private TreatmentMatcher() {
    }

    @NotNull
    public static MatchResult<BiopsyTreatmentData> matchTreatmentsToBiopsies(@NotNull final String patientId,
            @NotNull final List<BiopsyData> biopsies, @NotNull final List<BiopsyTreatmentData> treatments) {
        final List<BiopsyTreatmentData> matchedTreatments = Lists.newArrayList();
        final List<ValidationFinding> findings = Lists.newArrayList();

        List<BiopsyData> remainingBiopsies = biopsies;
        for (final BiopsyTreatmentData treatment : treatments) {
            final String treatmentGiven = treatment.treatmentGiven();
            if (treatmentGiven == null) {
                matchedTreatments.add(treatment);
                continue;
            }
            if (treatmentGiven.toLowerCase().equals("yes")) {
                final Map<Boolean, List<BiopsyData>> partitionsNoTreatment = remainingBiopsies.stream()
                        .collect(Collectors.partitioningBy(clinicalBiopsy -> isPossibleMatch(clinicalBiopsy.date(), treatment.startDate())));
                final List<BiopsyData> possibleMatchesWithTreatment = partitionsNoTreatment.get(true);
                if (possibleMatchesWithTreatment.size() == 0) {
                    findings.add(treatmentMatchFinding(patientId, "no biopsy match for treatment with treatment given is yes", treatment.toString()));
                    matchedTreatments.add(treatment);
                } else if (possibleMatchesWithTreatment.size() > 1) {
                    findings.add(treatmentMatchFinding(patientId,
                            "multiple biopsy matches for treatment with treatment given is yes",
                            treatment + ". biopsies:  " + possibleMatchesWithTreatment.stream().map(BiopsyData::toString).collect(Collectors.toList())));
                    matchedTreatments.add(treatment);
                } else if (possibleMatchesWithTreatment.get(0).date() == null) {
                    findings.add(treatmentMatchFinding(patientId,
                            "treatment matched biopsy with null date with treatment given is yes.",
                            treatment.toString()));
                    matchedTreatments.add(treatment);
                } else {
                    final BiopsyData clinicalBiopsy = possibleMatchesWithTreatment.get(0);
                    matchedTreatments.add(ImmutableBiopsyTreatmentData.builder().from(treatment).biopsyId(clinicalBiopsy.id()).build());
                    remainingBiopsies = partitionsNoTreatment.get(false);
                }
            } else if (treatmentGiven.toLowerCase().equals("no")) {
                final Map<Boolean, List<BiopsyData>> partitionsWithTreatment = remainingBiopsies.stream()
                        .collect(Collectors.partitioningBy(clinicalBiopsy -> isPossibleMatchNoTreatment(clinicalBiopsy.date(),
                                treatment.treatmentGiven())));
                final List<BiopsyData> possibleMatchesNoTreatment = partitionsWithTreatment.get(true);
                if (possibleMatchesNoTreatment.size() == 0) {
                    findings.add(treatmentMatchFinding(patientId,
                            "no biopsy match for treatment with treatment given is no",
                            treatment.toString()));
                    matchedTreatments.add(treatment);
                } else {
                    final BiopsyData clinicalBiopsy = possibleMatchesNoTreatment.get(0);
                    matchedTreatments.add(ImmutableBiopsyTreatmentData.builder().from(treatment).biopsyId(clinicalBiopsy.id()).build());
                    remainingBiopsies = partitionsWithTreatment.get(false);
                }
            }
        }
        return new MatchResult<>(matchedTreatments, findings);
    }

    private static boolean isPossibleMatch(@Nullable final LocalDate biopsyDate, @Nullable final LocalDate treatmentStartDate) {
        return biopsyDate == null || isWithinThreshold(biopsyDate, treatmentStartDate);
    }

    private static boolean isPossibleMatchNoTreatment(@Nullable final LocalDate biopsyDate, @Nullable final String treatmentGiven) {
        return biopsyDate == null || MatchingNoTreatment(biopsyDate, treatmentGiven);
    }

    private static boolean MatchingNoTreatment(@NotNull final LocalDate biopsyDate, @Nullable final String treatmentGiven) {
        return biopsyDate != null & treatmentGiven.toLowerCase().equals("no");
    }

    private static boolean isWithinThreshold(@NotNull final LocalDate biopsyDate, @Nullable final LocalDate treatmentStartDate) {
        return treatmentStartDate != null && (treatmentStartDate.isAfter(biopsyDate) || treatmentStartDate.isEqual(biopsyDate))
                && Duration.between(biopsyDate.atStartOfDay(), treatmentStartDate.atStartOfDay()).toDays()
                < Config.MAX_DAYS_BETWEEN_TREATMENT_AND_BIOPSY;
    }

    @NotNull
    private static ValidationFinding treatmentMatchFinding(@NotNull String patientIdentifier, @NotNull String message,
            @NotNull String details) {
        return ValidationFinding.of("match", patientIdentifier, FORM_TREATMENT, message, FormStatus.unknown(), details);
    }
}
