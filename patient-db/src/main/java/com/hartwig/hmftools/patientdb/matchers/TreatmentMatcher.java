package com.hartwig.hmftools.patientdb.matchers;

import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentReader.FORM_TREATMENT;

import java.time.Duration;
import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.Config;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class TreatmentMatcher {
    private static final Logger LOGGER = LogManager.getLogger(TreatmentMatcher.class);

    private TreatmentMatcher() {
    }

    @NotNull
    public static MatchResult<BiopsyTreatmentData> matchTreatmentsToBiopsies(@NotNull final String patientId,
            @NotNull final List<BiopsyData> biopsies, @NotNull final List<BiopsyTreatmentData> treatments) {
        final List<BiopsyTreatmentData> matchedTreatments = Lists.newArrayList();
        final List<ValidationFinding> findings = Lists.newArrayList();
        if (biopsies.size() < treatments.size()) {
            LOGGER.warn(patientId + ": has more treatments(" + treatments.size() + ") than biopsies(" + biopsies.size() + ").");
        }
        List<BiopsyData> remainingBiopsies = biopsies;
        for (final BiopsyTreatmentData treatment : treatments) {
            final String treatmentGiven = treatment.treatmentGiven();
            if (treatmentGiven == null || treatmentGiven.toLowerCase().equals("no")) {
                matchedTreatments.add(treatment);
                continue;
            }
            final Map<Boolean, List<BiopsyData>> partitions = remainingBiopsies.stream()
                    .collect(Collectors.partitioningBy(clinicalBiopsy -> isPossibleMatch(clinicalBiopsy.date(), treatment.startDate())));
            final List<BiopsyData> possibleMatches = partitions.get(true);
            if (possibleMatches.size() == 0) {
                findings.add(ValidationFinding.of("match", patientId, FORM_TREATMENT, "no biopsy match for treatment", "", "",
                        treatment.toString()));
                matchedTreatments.add(treatment);
            } else if (possibleMatches.size() > 1) {
                findings.add(ValidationFinding.of("match", patientId, FORM_TREATMENT, "multiple biopsy matches for treatment", "", "",
                        treatment + ". biopsies:  " + possibleMatches.stream().map(BiopsyData::toString).collect(Collectors.toList())));
                matchedTreatments.add(treatment);
            } else if (possibleMatches.size() == 1 && possibleMatches.get(0).date() == null) {
                findings.add(ValidationFinding.of("match", patientId, FORM_TREATMENT, "treatment matched biopsy with null date.", "", "",
                        treatment.toString()));
                matchedTreatments.add(treatment);
            } else {
                final BiopsyData clinicalBiopsy = possibleMatches.get(0);
                matchedTreatments.add(ImmutableBiopsyTreatmentData.of(treatment.id(), treatment.treatmentGiven(), treatment.startDate(),
                        treatment.endDate(), treatment.drugs(), clinicalBiopsy.id(), treatment.formStatus(), treatment.formLocked()));
                remainingBiopsies = partitions.get(false);
            }
        }
        return new MatchResult<>(matchedTreatments, findings);
    }

    private static boolean isPossibleMatch(@Nullable final LocalDate biopsyDate, @Nullable final LocalDate treatmentStartDate) {
        return biopsyDate == null || isWithinThreshold(biopsyDate, treatmentStartDate);
    }

    private static boolean isWithinThreshold(@Nullable final LocalDate biopsyDate, @Nullable final LocalDate treatmentStartDate) {
        return biopsyDate != null && treatmentStartDate != null && (treatmentStartDate.isAfter(biopsyDate) || treatmentStartDate.isEqual(
                biopsyDate)) && Duration.between(biopsyDate.atStartOfDay(), treatmentStartDate.atStartOfDay()).toDays()
                < Config.MAX_DAYS_BETWEEN_TREATMENT_AND_BIOPSY;
    }
}
