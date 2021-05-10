package com.hartwig.hmftools.patientdb.clinical.matchers;

import java.time.Duration;
import java.time.LocalDate;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ImmutableValidationFinding;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.clinical.ClinicalConstants;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableBiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;

import org.jetbrains.annotations.NotNull;

public final class BiopsyMatcher {

    private BiopsyMatcher() {
    }

    @NotNull
    public static MatchResult<BiopsyData> matchBiopsiesToTumorSamples(@NotNull String patientIdentifier,
            @NotNull List<SampleData> sequencedBiopsies, @NotNull List<BiopsyData> clinicalBiopsies) {
        List<BiopsyData> matchedBiopsies = Lists.newArrayList();
        List<ValidationFinding> findings = Lists.newArrayList();

        List<BiopsyData> remainingBiopsies = clinicalBiopsies;
        Collections.sort(remainingBiopsies);
        Collections.sort(sequencedBiopsies);

        if (clinicalBiopsies.size() < sequencedBiopsies.size()) {
            findings.add(biopsyMatchFinding(patientIdentifier,
                    "Not enough clinical biopsy forms to match for every sequenced sample",
                    "Clinical biopsies: " + clinicalBiopsies.size() + ", sequenced samples: " + sequencedBiopsies.size()));
        }

        for (SampleData sequencedBiopsy : sequencedBiopsies) {
            Map<Boolean, List<BiopsyData>> partitions = remainingBiopsies.stream()
                    .collect(Collectors.partitioningBy(clinicalBiopsy -> isPossibleMatch(sequencedBiopsy, clinicalBiopsy)));
            List<BiopsyData> possibleMatches = partitions.get(true);
            if (possibleMatches.size() == 1 && possibleMatches.get(0).date() != null) {
                BiopsyData clinicalBiopsy = possibleMatches.get(0);
                matchedBiopsies.add(ImmutableBiopsyData.builder().from(clinicalBiopsy).sampleId(sequencedBiopsy.sampleId()).build());
                if (hasMismatchOnSamplingDate(sequencedBiopsy, clinicalBiopsy) && !patientIdentifier.startsWith("WIDE")) {
                    findings.add(biopsyMatchFinding(patientIdentifier,
                            "Sampling date does not equal biopsy date in matched biopsy",
                            "sampling date: " + sequencedBiopsy.samplingDate() + "; ecrf biopsy date: " + clinicalBiopsy.date()));
                }
                remainingBiopsies = partitions.get(false);
            } else {
                if (possibleMatches.size() == 0 || (possibleMatches.size() == 1 && possibleMatches.get(0).date() == null)) {
                    findings.add(biopsyMatchFinding(patientIdentifier,
                            "Could not match any clinical biopsy with sequenced sample.",
                            "sample: " + sampleDataToString(sequencedBiopsy) + "; ecrf biopsies: " + clinicalBiopsies + ". match criteria: "
                                    + getMatchDateCriteria(sequencedBiopsy)));
                } else if (possibleMatches.size() > 1) {
                    findings.add(biopsyMatchFinding(patientIdentifier,
                            "More than 1 possible clinical biopsy match for sequenced sample.",
                            "sample: " + sampleDataToString(sequencedBiopsy) + "; ecrf biopsies: " + clinicalBiopsies + ". match criteria: "
                                    + getMatchDateCriteria(sequencedBiopsy)));
                } else {
                    findings.add(biopsyMatchFinding(patientIdentifier,
                            "Undetermined match issue in biopsy matcher",
                            "sample: " + sampleDataToString(sequencedBiopsy)));
                }
            }
        }
        matchedBiopsies.addAll(remainingBiopsies);
        return new MatchResult<>(matchedBiopsies, findings);
    }

    @NotNull
    private static String sampleDataToString(@NotNull SampleData sampleData) {
        return sampleData.sampleId() + "(S:" + sampleData.samplingDate() + ", A:" + sampleData.arrivalDate() + ")";
    }

    private static boolean isPossibleMatch(@NotNull SampleData sequencedBiopsy, @NotNull BiopsyData clinicalBiopsy) {
        return clinicalBiopsy.date() == null || (isWithinThreshold(sequencedBiopsy, clinicalBiopsy)
                && clinicalBiopsy.isPotentiallyEvaluable());
    }

    private static boolean isWithinThreshold(@NotNull SampleData sequencedBiopsy, @NotNull BiopsyData clinicalBiopsy) {
        LocalDate biopsyDate = clinicalBiopsy.date();
        assert biopsyDate != null;

        LocalDate samplingDate = sequencedBiopsy.samplingDate();
        if (samplingDate != null) {
            return Math.abs(Duration.between(biopsyDate.atStartOfDay(), samplingDate.atStartOfDay()).toDays())
                    < ClinicalConstants.MAX_DAYS_BETWEEN_SAMPLING_AND_BIOPSY_DATE;
        } else {
            long daysFromBiopsyTakenToArrival =
                    Duration.between(biopsyDate.atStartOfDay(), sequencedBiopsy.arrivalDate().atStartOfDay()).toDays();
            return daysFromBiopsyTakenToArrival >= 0
                    && daysFromBiopsyTakenToArrival < ClinicalConstants.MAX_DAYS_ARRIVAL_DATE_AFTER_BIOPSY_DATE;
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
    private static String getMatchDateCriteria(@NotNull SampleData sampleData) {
        if (sampleData.samplingDate() != null) {
            return "Sampling date " + sampleData.samplingDate() + " threshold: "
                    + ClinicalConstants.MAX_DAYS_BETWEEN_SAMPLING_AND_BIOPSY_DATE;
        }
        return "Arrival date " + sampleData.arrivalDate() + " threshold: " + ClinicalConstants.MAX_DAYS_ARRIVAL_DATE_AFTER_BIOPSY_DATE;
    }

    @NotNull
    private static ValidationFinding biopsyMatchFinding(@NotNull String patientIdentifier, @NotNull String finding,
            @NotNull String details) {
        return ImmutableValidationFinding.builder()
                .level("match")
                .patientIdentifier(patientIdentifier)
                .message(finding)
                .formStatus(FormStatus.undefined())
                .details(details)
                .build();
    }
}
