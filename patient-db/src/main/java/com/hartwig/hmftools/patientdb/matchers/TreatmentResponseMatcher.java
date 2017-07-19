package com.hartwig.hmftools.patientdb.matchers;

import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FORM_TUMOR_MEASUREMENT;

import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ImmutableValidationFinding;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentResponseData;

import org.jetbrains.annotations.NotNull;

public final class TreatmentResponseMatcher {
    private TreatmentResponseMatcher() {
    }

    public static MatchResult<BiopsyTreatmentResponseData> matchTreatmentResponsesToTreatments(@NotNull final String patientId,
            @NotNull final List<BiopsyTreatmentData> treatments, @NotNull final List<BiopsyTreatmentResponseData> responses) {
        final List<BiopsyTreatmentResponseData> matchedResponses = Lists.newArrayList();
        final List<ValidationFinding> findings = Lists.newArrayList();
        for (final BiopsyTreatmentResponseData response : responses) {
            final Map<Boolean, List<BiopsyTreatmentData>> partitions =
                    treatments.stream().collect(Collectors.partitioningBy(treatment -> responseMatchesTreatment(response, treatment)));
            final List<BiopsyTreatmentData> matchedTreatments = partitions.get(true);
            if (matchedTreatments.size() == 0) {
                matchedResponses.add(response);
            } else if (matchedTreatments.size() > 1) {
                findings.add(new ImmutableValidationFinding("match", patientId, FORM_TUMOR_MEASUREMENT,
                        "treatment response " + response.response() + "(" + response.date() + ") matches multiple treatments:"
                                + matchedTreatments.stream()
                                .map(treatment -> treatment.treatmentName() + "(" + treatment.startDate() + " - " + treatment.endDate()
                                        + ")")
                                .collect(Collectors.toList()), "", ""));
                matchedResponses.add(response);
            } else {
                matchedResponses.add(ImmutableBiopsyTreatmentResponseData.of(matchedTreatments.get(0).id(), response.assessmentDate(),
                        response.responseDate(), response.response(), response.measurementDone(), response.formStatus(),
                        response.formLocked()));
            }
        }
        return new MatchResult<>(matchedResponses, findings);
    }

    public static boolean responseMatchesTreatment(@NotNull final BiopsyTreatmentResponseData response,
            @NotNull final BiopsyTreatmentData treatment) {
        final LocalDate treatmentStart = treatment.startDate();
        final LocalDate treatmentEnd = treatment.endDate();
        final LocalDate responseDate = response.date();
        return !(treatmentStart == null || responseDate == null) && (responseDate.isAfter(treatmentStart) && (treatmentEnd == null
                || responseDate.isBefore(treatmentEnd) || responseDate.isEqual(treatmentEnd)));
    }
}
