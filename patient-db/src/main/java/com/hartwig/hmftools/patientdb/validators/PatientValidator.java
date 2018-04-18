package com.hartwig.hmftools.patientdb.validators;

import static java.util.Comparator.comparing;
import static java.util.Comparator.naturalOrder;
import static java.util.Comparator.nullsLast;

import java.time.Duration;
import java.time.LocalDate;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.Config;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.CuratedTreatment;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.TreatmentData;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PatientValidator {

    private static final String ECRF_LEVEL = "ecrf";

    @NotNull
    public static List<ValidationFinding> validatePatient(@NotNull final Patient patient) {
        final String patientIdentifier = patient.patientIdentifier();

        final List<ValidationFinding> findings = Lists.newArrayList();

        findings.addAll(validateBaselineData(patientIdentifier, patient.baselineData()));
        findings.addAll(validatePreTreatmentData(patientIdentifier, patient.preTreatmentData()));
        findings.addAll(validateBiopsies(patientIdentifier, patient.clinicalBiopsies(), patient.sequencedBiopsies().size()));
        findings.addAll(validateTreatments(patientIdentifier, patient.treatments()));
        findings.addAll(validateTreatmentResponses(patientIdentifier, patient.treatments(), patient.treatmentResponses()));
        findings.addAll(validateDeathDate(patientIdentifier, patient.baselineData(), patient.treatments()));
        findings.addAll(validateInformedConsentDate(patientIdentifier, patient.baselineData(), patient.clinicalBiopsies()));

        findings.addAll(validateTumorLocationCuration(patientIdentifier, patient.baselineData()));
        findings.addAll(validateBiopsyTreatmentsCuration(patientIdentifier, patient.treatments()));
        findings.addAll(validatePreTreatmentCuration(patientIdentifier, patient.preTreatmentData()));

        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateBaselineData(@NotNull String patientIdentifier, @NotNull BaselineData baselineData) {
        final List<ValidationFinding> findings = Lists.newArrayList();
        if (baselineData.curatedTumorLocation().searchTerm() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "primary tumor location empty",
                    baselineData.primaryTumorStatus()));
        }

        if (baselineData.gender() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "gender empty", baselineData.demographyStatus()));
        }

        if (baselineData.registrationDate() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "registration date empty or in wrong format",
                    FormStatus.merge(baselineData.selectionCriteriaStatus(), baselineData.eligibilityStatus())));
        }

        if (baselineData.informedConsentDate() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "informed consent date empty or in wrong format",
                    baselineData.informedConsentStatus()));
        }

        if (baselineData.birthYear() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "birth year could not be determined",
                    FormStatus.merge(baselineData.eligibilityStatus(), baselineData.selectionCriteriaStatus())));
        }

        if (baselineData.hospital() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "hospital could not be determined",
                    FormStatus.merge(baselineData.eligibilityStatus(), baselineData.selectionCriteriaStatus())));
        }

        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validatePreTreatmentData(@NotNull final String patientIdentifier,
            @NotNull final PreTreatmentData preTreatmentData) {
        final String preTreatmentGiven = preTreatmentData.treatmentGiven();
        final String preRadioTreatmentGiven = preTreatmentData.radiotherapyGiven();
        final List<ValidationFinding> findings = Lists.newArrayList();

        if (preTreatmentGiven == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "pre systemic treatment given empty",
                    preTreatmentData.formStatus()));
        }

        if (preRadioTreatmentGiven == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "pre radio treatment given empty",
                    preTreatmentData.formStatus()));
        }

        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateBiopsies(@NotNull final String patientIdentifier, @NotNull final List<BiopsyData> biopsies,
            int expectedMatches) {
        final List<ValidationFinding> findings = Lists.newArrayList();

        biopsies.forEach(biopsy -> findings.addAll(validateBiopsyData(patientIdentifier, biopsy)));

        if (biopsies.size() < expectedMatches) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "not enough biopsy forms to match with sequenced samples",
                    FormStatus.undefined()));
        }

        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateBiopsyData(@NotNull final String patientIdentifier, @NotNull final BiopsyData biopsyData) {
        final List<ValidationFinding> findings = Lists.newArrayList();

        if (biopsyData.date() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "biopsy date empty or in wrong format",
                    biopsyData.formStatus()));
        }

        if (biopsyData.site() == null && biopsyData.location() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "biopsy site and biopsy location are empty",
                    biopsyData.formStatus()));
        }

        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateTreatments(@NotNull final String patientIdentifier,
            @NotNull final List<BiopsyTreatmentData> treatments) {
        final List<ValidationFinding> findings = Lists.newArrayList();
        treatments.forEach(treatment -> findings.addAll(validateTreatmentData(patientIdentifier, treatment)));
        Collections.sort(treatments);

        if (treatments.size() > 1) {
            for (int index = 1; index < treatments.size(); index++) {
                final BiopsyTreatmentData currentTreatment = treatments.get(index);
                final LocalDate startDate = currentTreatment.startDate();
                final LocalDate previousTreatmentEnd = treatments.get(index - 1).endDate();
                if (startDate != null && (previousTreatmentEnd == null || startDate.isBefore(previousTreatmentEnd))) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "subsequent treatment starts before the end of previous treatment",
                            currentTreatment.formStatus()));
                }
            }
            final List<BiopsyTreatmentData> nonFinishedTreatments =
                    treatments.stream().filter(treatment -> treatment.endDate() == null).collect(Collectors.toList());
            if (nonFinishedTreatments.size() > 1) {
                nonFinishedTreatments.forEach(treatment -> findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "end of at least 1 non-final treatment is missing",
                        treatment.formStatus())));
            }
        }
        return findings;
    }

    @NotNull
    private static List<ValidationFinding> validateTreatmentData(@NotNull final String patientIdentifier,
            @NotNull final BiopsyTreatmentData treatmentData) {
        final String treatmentGiven = treatmentData.treatmentGiven();
        final List<ValidationFinding> findings = Lists.newArrayList();
        if (treatmentGiven == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "treatment given field empty", treatmentData.formStatus()));
        } else if (treatmentGiven.trim().toLowerCase().equals("yes")) {
            if (treatmentData.drugs().isEmpty()) {
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "treatment given is yes, but no treatment data filled in",
                        treatmentData.formStatus()));
            } else {
                treatmentData.drugs()
                        .forEach(drug -> findings.addAll(validateDrugData(patientIdentifier, drug, treatmentData.formStatus())));
            }
        } else if (treatmentGiven.trim().toLowerCase().equals("no")) {
            if (!treatmentData.drugs().isEmpty()) {
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "treatment given is no, but treatment data is filled in",
                        treatmentData.formStatus()));
            }
        } else {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "treatment given is not yes/no", treatmentData.formStatus()));
        }

        final String postRadioTherapy = treatmentData.radiotherapyGiven();
        if (postRadioTherapy == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "radio therapy given field empty",
                    treatmentData.formStatus()));
        }
        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateDrugData(@NotNull final String patientIdentifier, @NotNull final DrugData drugData,
            @NotNull final FormStatus formStatus) {
        final LocalDate drugStart = drugData.startDate();
        final List<ValidationFinding> findings = Lists.newArrayList();
        if (drugStart == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "drug start date empty or in wrong format", formStatus));
        } else {
            final LocalDate drugEnd = drugData.endDate();
            if (drugEnd != null) {
                if (drugStart.isAfter(drugEnd)) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "drug startDate is after drug endDate", formStatus));
                }
            }
        }
        if (drugData.name() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "drug name empty", formStatus));
        }
        return findings;
    }

    @NotNull
    private static List<ValidationFinding> validatePreTreatmentCuration(@NotNull final String patientIdentifier,
            @NotNull PreTreatmentData preTreatmentData) {
        return validateTreatmentCuration(patientIdentifier, "preTreatmentCuration", Lists.newArrayList(preTreatmentData));
    }

    @NotNull
    private static List<ValidationFinding> validateBiopsyTreatmentsCuration(@NotNull final String patientIdentifier,
            @NotNull List<BiopsyTreatmentData> treatments) {
        return validateTreatmentCuration(patientIdentifier, "treatmentCuration", treatments);
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateTreatmentCuration(@NotNull final String patientIdentifier, @NotNull String curationName,
            @NotNull final List<? extends TreatmentData> treatments) {
        final List<ValidationFinding> findings = Lists.newArrayList();
        treatments.forEach(treatment -> treatment.drugs().forEach(drug -> {
            final String drugName = drug.name();
            if (drugName != null) {
                if (drug.curatedTreatments().isEmpty()) {
                    findings.add(ValidationFinding.of(curationName,
                            patientIdentifier,
                            "failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous",
                            treatment.formStatus(),
                            drugName));
                } else {
                    final List<String> curatedTreatments =
                            drug.curatedTreatments().stream().map(CuratedTreatment::name).collect(Collectors.toList());
                    final List<String> matchedTerms =
                            drug.curatedTreatments().stream().map(CuratedTreatment::searchTerm).collect(Collectors.toList());
                    final long lengthOfMatchedCharacters = matchedTerms.stream().mapToLong(String::length).sum();
                    final long lengthOfSearchCharacters = drugName.chars().filter(Character::isLetterOrDigit).count();
                    if (lengthOfMatchedCharacters > 0 && (double) lengthOfMatchedCharacters / lengthOfSearchCharacters < 0.9) {
                        findings.add(ValidationFinding.of(curationName,
                                patientIdentifier,
                                "matched drugs are based on less than 90% of search term.",
                                treatment.formStatus(),
                                drugName + " matched to " + Strings.join(curatedTreatments, ',') + " based on " + Strings.join(matchedTerms,
                                        ',')));
                    }
                }
            }
        }));
        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateTumorLocationCuration(@NotNull String patientIdentifier, @NotNull BaselineData baselineData) {
        final List<ValidationFinding> findings = Lists.newArrayList();
        final String searchTerm = baselineData.curatedTumorLocation().searchTerm();
        if (searchTerm != null && baselineData.curatedTumorLocation().primaryTumorLocation() == null) {
            findings.add(ValidationFinding.of("tumorLocationCuration",
                    patientIdentifier,
                    "failed to curate primary tumor location",
                    baselineData.primaryTumorStatus(),
                    searchTerm));
        }
        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateTreatmentResponses(@NotNull final String patientIdentifier,
            @NotNull final List<BiopsyTreatmentData> treatments, @NotNull final List<BiopsyTreatmentResponseData> responses) {
        final List<ValidationFinding> findings = Lists.newArrayList();
        for (int i = 0; i < responses.size(); i++) {
            findings.addAll(validateTreatmentResponse(patientIdentifier, responses.get(i), i == 0));
        }

        Collections.sort(treatments);
        Collections.sort(responses);
        if (treatments.isEmpty() && !responses.isEmpty()) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "treatment response filled in, but treatment data missing",
                    FormStatus.undefined()));
        }
        if (!treatments.isEmpty() && !responses.isEmpty()) {
            final LocalDate firstResponseDate = responses.get(0).date();
            final LocalDate firstTreatmentStart = treatments.get(0).startDate();
            if (firstResponseDate != null && firstTreatmentStart != null && firstResponseDate.isAfter(firstTreatmentStart)) {
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "first (baseline) measurement date is after first treatment start",
                        FormStatus.merge(treatments.get(0).formStatus(), responses.get(0).formStatus()),
                        "first treatment response: " + firstResponseDate + "; first treatment start: " + firstTreatmentStart));

            }
        }
        final List<BiopsyTreatmentData> treatmentsMissingResponse = treatments.stream()
                .filter(treatment -> shouldHaveResponse(treatment) && !hasResponse(treatment.id(), responses))
                .collect(Collectors.toList());

        if (treatmentsMissingResponse.size() > 0) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "no treatment response for at least 1 treatment",
                    FormStatus.undefined(),
                    "treatments " + treatmentsMissingResponse.stream().map(BiopsyTreatmentData::toString).collect(Collectors.toList())
                            + " should have response since they lasted more than " + Config.IMMEDIATE_TREATMENT_END_THRESHOLD
                            + " days and started more than " + Config.START_DATE_RESPONSE_THRESHOLD + " days ago"));
        }
        return findings;
    }

    private static boolean shouldHaveResponse(@NotNull final BiopsyTreatmentData treatment) {
        final LocalDate treatmentStart = treatment.startDate();
        return (treatmentStart != null && Duration.between(treatmentStart.atStartOfDay(), LocalDate.now().atStartOfDay()).toDays()
                > Config.START_DATE_RESPONSE_THRESHOLD) && !endedImmediately(treatment);
    }

    private static boolean endedImmediately(@NotNull final BiopsyTreatmentData treatment) {
        final LocalDate treatmentStart = treatment.startDate();
        final LocalDate treatmentEnd = treatment.endDate();
        return treatmentStart != null && treatmentEnd != null
                && Duration.between(treatmentStart.atStartOfDay(), treatmentEnd.atStartOfDay()).toDays()
                < Config.IMMEDIATE_TREATMENT_END_THRESHOLD;
    }

    private static boolean hasResponse(@NotNull Integer treatmentId, @NotNull List<BiopsyTreatmentResponseData> responses) {
        return responses.stream().anyMatch(response -> treatmentId.equals(response.treatmentId()));
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateTreatmentResponse(@NotNull final String patientIdentifier,
            @NotNull final BiopsyTreatmentResponseData treatmentResponse, boolean isFirstResponse) {
        final List<ValidationFinding> findings = Lists.newArrayList();
        final String measurementDone = treatmentResponse.measurementDone();
        final String response = treatmentResponse.response();
        final LocalDate date = treatmentResponse.date();

        if (measurementDone == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "measurement done field empty",
                    treatmentResponse.formStatus()));

        } else if (measurementDone.trim().toLowerCase().equals("yes")) {
            if (date == null) {
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "response date and assessment date empty or in wrong format",
                        treatmentResponse.formStatus()));
            }
            if (response == null && !isFirstResponse) {
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "measurement done is yes, but response is empty (non-first response)",
                        treatmentResponse.formStatus()));
            }
        } else if (measurementDone.trim().equalsIgnoreCase("no")) {
            if (!treatmentResponse.isNotDoneResponse()) {
                if (date != null) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "measurement done is no, but assessment date or response date is filled in",
                            treatmentResponse.formStatus(),
                            "effective response date: " + date));
                }
                if (response != null) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "measurement done is no, but response filled in",
                            treatmentResponse.formStatus(),
                            "response: " + response));
                }
            }
        } else {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "measurement done is not yes/no",
                    treatmentResponse.formStatus(),
                    "measurement done: " + measurementDone));
        }
        if (response != null && date == null && !treatmentResponse.isNotDoneResponse()) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "response filled in, but no assessment date and response date found",
                    treatmentResponse.formStatus(),
                    "response: " + response));
        }
        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateDeathDate(@NotNull final String patientIdentifier, @NotNull final BaselineData baselineData,
            @NotNull final List<BiopsyTreatmentData> treatments) {
        final List<ValidationFinding> findings = Lists.newArrayList();
        final LocalDate deathDate = baselineData.deathDate();
        if (deathDate != null && !treatments.isEmpty()) {
            treatments.sort(comparing(BiopsyTreatmentData::endDate, nullsLast(naturalOrder())));
            final BiopsyTreatmentData lastTreatment = treatments.get(treatments.size() - 1);
            final String treatmentGiven = lastTreatment.treatmentGiven();
            if (treatmentGiven != null && treatmentGiven.equalsIgnoreCase("yes")) {
                final LocalDate lastTreatmentEndDate = lastTreatment.endDate();
                final LocalDate firstTreatmentStart = treatments.get(0).startDate();
                if (lastTreatmentEndDate == null || lastTreatmentEndDate.isAfter(deathDate)) {
                    String details = "death date (" + deathDate + ") before end of last treatment (" + lastTreatmentEndDate + ")"
                            + " and start first treatment is (" + firstTreatmentStart + ")" + " and treatmentGiven: " + treatmentGiven
                            + ")";

                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "death date before end of last treatment",
                            FormStatus.merge(baselineData.deathStatus(), lastTreatment.formStatus()),
                            details));
                }
            }

        }
        return findings;
    }

    @NotNull
    static List<ValidationFinding> validateInformedConsentDate(@NotNull final String patientIdentifier,
            @NotNull final BaselineData baselineData, @NotNull final List<BiopsyData> biopsies) {
        final List<ValidationFinding> findings = Lists.newArrayList();
        final LocalDate informedConsentDate = baselineData.informedConsentDate();
        if (informedConsentDate != null && !biopsies.isEmpty()) {
            final List<BiopsyData> biopsiesPriorToInformedConsent = biopsies.stream().filter(biopsy -> {
                final LocalDate biopsyDate = biopsy.date();
                return biopsyDate != null && biopsyDate.isBefore(informedConsentDate);
            }).collect(Collectors.toList());
            if (biopsiesPriorToInformedConsent.size() > 0) {
                final String detailsMessage =
                        "informedConsentDate: " + informedConsentDate + ". biopsies: " + biopsiesPriorToInformedConsent.stream()
                                .map(BiopsyData::toString)
                                .collect(Collectors.toList());

                FormStatus mergedFormStatus = baselineData.informedConsentStatus();
                for (BiopsyData biopsy : biopsies) {
                    mergedFormStatus = FormStatus.merge(mergedFormStatus, biopsy.formStatus());
                }
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "at least 1 biopsy taken before informed consent date",
                        mergedFormStatus,
                        detailsMessage));
            }
        }
        return findings;
    }
}
