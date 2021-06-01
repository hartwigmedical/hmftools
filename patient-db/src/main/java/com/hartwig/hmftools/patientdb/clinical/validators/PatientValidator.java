package com.hartwig.hmftools.patientdb.clinical.validators;

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
import com.hartwig.hmftools.patientdb.clinical.ClinicalConstants;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedDrug;
import com.hartwig.hmftools.patientdb.clinical.datamodel.DrugData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.PreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.TreatmentData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PatientValidator {

    private static final String ECRF_LEVEL = "ecrf";

    private PatientValidator() {
    }

    @NotNull
    public static List<ValidationFinding> validatePatient(@NotNull Patient patient) {
        String patientIdentifier = patient.patientIdentifier();

        List<ValidationFinding> findings = Lists.newArrayList();

        findings.addAll(validateBaselineData(patientIdentifier, patient.baselineData()));
        findings.addAll(validatePreTreatmentData(patientIdentifier, patient.preTreatmentData()));
        findings.addAll(validateBiopsies(patientIdentifier, patient.clinicalBiopsies()));
        findings.addAll(validateTreatments(patientIdentifier, patient.treatments()));
        findings.addAll(validateTreatmentResponses(patientIdentifier, patient.treatments(), patient.treatmentResponses()));
        findings.addAll(validateDeathDate(patientIdentifier, patient.baselineData(), patient.treatments()));
        findings.addAll(validateInformedConsentDate(patientIdentifier, patient.baselineData(), patient.clinicalBiopsies()));

        findings.addAll(validatePrimaryTumorCuration(patientIdentifier, patient.baselineData()));
        findings.addAll(validateBiopsyTreatmentsCuration(patientIdentifier, patient.treatments()));
        findings.addAll(validatePreTreatmentCuration(patientIdentifier, patient.preTreatmentData()));

        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateBaselineData(@NotNull String patientIdentifier, @NotNull BaselineData baselineData) {
        List<ValidationFinding> findings = Lists.newArrayList();
        if (baselineData.curatedPrimaryTumor().searchTerm() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Primary tumor location empty",
                    baselineData.primaryTumorStatus()));
        }

        if (baselineData.gender() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "Gender empty", baselineData.demographyStatus()));
        }

        if (baselineData.registrationDate() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Registration date empty or in wrong format",
                    FormStatus.merge(baselineData.selectionCriteriaStatus(), baselineData.eligibilityStatus())));
        }

        if (baselineData.informedConsentDate() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Informed consent date empty or in wrong format",
                    baselineData.informedConsentStatus()));
        }

        if (baselineData.birthYear() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Birth year could not be determined",
                    FormStatus.merge(baselineData.eligibilityStatus(), baselineData.selectionCriteriaStatus())));
        }

        if (baselineData.hospital() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Hospital could not be determined",
                    FormStatus.merge(baselineData.eligibilityStatus(), baselineData.selectionCriteriaStatus())));
        }

        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validatePreTreatmentData(@NotNull String patientIdentifier, @NotNull PreTreatmentData preTreatmentData) {
        String preTreatmentGiven = preTreatmentData.treatmentGiven();
        String preRadioTreatmentGiven = preTreatmentData.radiotherapyGiven();
        List<ValidationFinding> findings = Lists.newArrayList();

        if (preTreatmentGiven == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Pre systemic treatment given empty",
                    preTreatmentData.formStatus()));
        }

        if (preRadioTreatmentGiven == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Pre radio treatment given empty",
                    preTreatmentData.formStatus()));
        }

        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateBiopsies(@NotNull String patientIdentifier, @NotNull List<BiopsyData> biopsies) {
        List<ValidationFinding> findings = Lists.newArrayList();

        for (BiopsyData biopsy : biopsies) {
            if (biopsy.sampleId() != null) {
                if (biopsy.date() == null) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "Biopsy date empty or in wrong format",
                            biopsy.formStatus()));
                }

                if (biopsy.site() == null && biopsy.location() == null) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "Biopsy site and biopsy location are empty",
                            biopsy.formStatus()));
                }
            }
        }

        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateTreatments(@NotNull String patientIdentifier, @NotNull List<BiopsyTreatmentData> treatments) {
        List<ValidationFinding> findings = Lists.newArrayList();
        List<BiopsyTreatmentData> matchedTreatments =
                treatments.stream().filter(treatment -> treatment.biopsyId() != null).collect(Collectors.toList());
        matchedTreatments.forEach(treatment -> findings.addAll(validateTreatmentData(patientIdentifier, treatment)));
        Collections.sort(matchedTreatments);

        if (matchedTreatments.size() > 1) {
            for (int index = 1; index < matchedTreatments.size(); index++) {
                BiopsyTreatmentData currentTreatment = matchedTreatments.get(index);
                LocalDate startDate = currentTreatment.startDate();
                LocalDate previousTreatmentEnd = matchedTreatments.get(index - 1).endDate();
                if (startDate != null && (previousTreatmentEnd == null || startDate.isBefore(previousTreatmentEnd))) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "Subsequent treatment starts before the end of previous treatment",
                            currentTreatment.formStatus()));
                }
            }
            List<BiopsyTreatmentData> nonFinishedTreatments = matchedTreatments.stream()
                    .filter(treatment -> treatment.startDate() != null)
                    .filter(treatment -> treatment.endDate() == null)
                    .collect(Collectors.toList());
            if (nonFinishedTreatments.size() > 1) {
                nonFinishedTreatments.forEach(treatment -> findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "End of at least 1 non-final treatment is missing",
                        treatment.formStatus())));
            }
        }
        return findings;
    }

    @NotNull
    private static List<ValidationFinding> validateTreatmentData(@NotNull String patientIdentifier,
            @NotNull BiopsyTreatmentData treatment) {
        String treatmentGiven = treatment.treatmentGiven();
        List<ValidationFinding> findings = Lists.newArrayList();
        if (treatmentGiven == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "Treatment given field empty", treatment.formStatus()));
        } else if (treatmentGiven.equalsIgnoreCase("yes")) {
            if (treatment.drugs().isEmpty()) {
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "Treatment given is yes, but no drug are filled in",
                        treatment.formStatus()));
            } else {
                treatment.drugs().forEach(drug -> findings.addAll(validateDrugData(patientIdentifier, drug, treatment.formStatus())));
            }
        } else if (treatmentGiven.equalsIgnoreCase("no")) {
            if (!treatment.drugs().isEmpty()) {
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "Treatment given is no, but drug are filled in",
                        treatment.formStatus(),
                        treatment.drugs().toString()));
            }
        } else {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Treatment given is not yes/no",
                    treatment.formStatus(),
                    treatmentGiven));
        }

        String postRadioTherapy = treatment.radiotherapyGiven();
        if (postRadioTherapy == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "Radio therapy given field empty", treatment.formStatus()));
        }
        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateDrugData(@NotNull String patientIdentifier, @NotNull DrugData drugData,
            @NotNull FormStatus formStatus) {
        LocalDate drugStart = drugData.startDate();
        List<ValidationFinding> findings = Lists.newArrayList();
        if (drugStart == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "Drug start date empty or in wrong format", formStatus));
        } else {
            LocalDate drugEnd = drugData.endDate();
            if (drugEnd != null) {
                if (drugStart.isAfter(drugEnd)) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "Drug start date is after drug end date", formStatus));
                }
            }
        }
        if (drugData.name() == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL, patientIdentifier, "Drug name empty", formStatus));
        }
        return findings;
    }

    @NotNull
    private static List<ValidationFinding> validatePreTreatmentCuration(@NotNull String patientIdentifier,
            @NotNull PreTreatmentData preTreatmentData) {
        return validateTreatmentCuration(patientIdentifier, "preTreatmentCuration", Lists.newArrayList(preTreatmentData));
    }

    @NotNull
    private static List<ValidationFinding> validateBiopsyTreatmentsCuration(@NotNull String patientIdentifier,
            @NotNull List<BiopsyTreatmentData> treatments) {
        return validateTreatmentCuration(patientIdentifier, "treatmentCuration", treatments);
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateTreatmentCuration(@NotNull String patientIdentifier, @NotNull String curationName,
            @NotNull List<? extends TreatmentData> treatments) {
        List<ValidationFinding> findings = Lists.newArrayList();
        treatments.forEach(treatment -> treatment.drugs().forEach(drug -> {
            String drugName = drug.name();
            if (drugName != null) {
                if (drug.curatedDrugs().isEmpty()) {
                    findings.add(ValidationFinding.of(curationName,
                            patientIdentifier,
                            "Failed to curate ecrf drug. Curated list contained no matching entry, or match was ambiguous",
                            treatment.formStatus(),
                            drugName));
                } else {
                    List<String> curatedTreatments = drug.curatedDrugs().stream().map(CuratedDrug::name).collect(Collectors.toList());
                    List<String> matchedTerms = drug.curatedDrugs().stream().map(CuratedDrug::searchTerm).collect(Collectors.toList());
                    long lengthOfMatchedCharacters = matchedTerms.stream().mapToLong(String::length).sum();
                    long lengthOfSearchCharacters = drugName.chars().filter(Character::isLetterOrDigit).count();
                    if (lengthOfMatchedCharacters > 0 && (double) lengthOfMatchedCharacters / lengthOfSearchCharacters < 0.9) {
                        findings.add(ValidationFinding.of(curationName,
                                patientIdentifier,
                                "Matched drugs are based on less than 90% of search term.",
                                treatment.formStatus(),
                                "'" + drugName + "' matched to '" + Strings.join(curatedTreatments, ',') + "' based on '" + Strings.join(
                                        matchedTerms,
                                        ',') + "'"));
                    }
                }
            }
        }));
        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validatePrimaryTumorCuration(@NotNull String patientIdentifier, @NotNull BaselineData baselineData) {
        List<ValidationFinding> findings = Lists.newArrayList();

        String searchTerm = baselineData.curatedPrimaryTumor().searchTerm();
        if (searchTerm != null && !searchTerm.isEmpty() && baselineData.curatedPrimaryTumor().location() == null) {
            findings.add(ValidationFinding.of("primaryTumorCuration",
                    patientIdentifier,
                    "Failed to curate primary tumor",
                    baselineData.primaryTumorStatus(),
                    searchTerm));
        }
        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateTreatmentResponses(@NotNull String patientIdentifier,
            @NotNull List<BiopsyTreatmentData> treatments, @NotNull List<BiopsyTreatmentResponseData> responses) {
        List<ValidationFinding> findings = Lists.newArrayList();
        for (int i = 0; i < responses.size(); i++) {
            findings.addAll(validateTreatmentResponse(patientIdentifier, responses.get(i), i == 0));
        }

        List<BiopsyTreatmentData> matchedTreatmentsMissingResponse = treatments.stream()
                .filter(treatment -> shouldHaveResponse(treatment) && !hasResponse(treatment.id(), responses)
                        && treatment.biopsyId() != null)
                .collect(Collectors.toList());

        if (matchedTreatmentsMissingResponse.size() > 0) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "No treatment response for at least 1 matched treatment",
                    FormStatus.undefined(),
                    "treatments " + matchedTreatmentsMissingResponse.stream()
                            .map(BiopsyTreatmentData::toString)
                            .collect(Collectors.toList()) + " should have response since they lasted more than "
                            + ClinicalConstants.IMMEDIATE_TREATMENT_END_THRESHOLD + " days and started more than "
                            + ClinicalConstants.START_DATE_RESPONSE_THRESHOLD + " days ago"));
        }
        return findings;
    }

    private static boolean shouldHaveResponse(@NotNull BiopsyTreatmentData treatment) {
        LocalDate treatmentStart = treatment.startDate();
        return (treatmentStart != null && Duration.between(treatmentStart.atStartOfDay(), LocalDate.now().atStartOfDay()).toDays()
                > ClinicalConstants.START_DATE_RESPONSE_THRESHOLD) && !endedImmediately(treatment);
    }

    private static boolean endedImmediately(@NotNull BiopsyTreatmentData treatment) {
        LocalDate treatmentStart = treatment.startDate();
        LocalDate treatmentEnd = treatment.endDate();
        return treatmentStart != null && treatmentEnd != null
                && Duration.between(treatmentStart.atStartOfDay(), treatmentEnd.atStartOfDay()).toDays()
                < ClinicalConstants.IMMEDIATE_TREATMENT_END_THRESHOLD;
    }

    private static boolean hasResponse(@NotNull Integer treatmentId, @NotNull List<BiopsyTreatmentResponseData> responses) {
        return responses.stream().anyMatch(response -> treatmentId.equals(response.treatmentId()));
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateTreatmentResponse(@NotNull String patientIdentifier,
            @NotNull BiopsyTreatmentResponseData treatmentResponse, boolean isFirstResponse) {
        List<ValidationFinding> findings = Lists.newArrayList();
        String measurementDone = treatmentResponse.measurementDone();
        String response = treatmentResponse.response();
        LocalDate date = treatmentResponse.date();

        if (measurementDone == null) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Measurement done field empty",
                    treatmentResponse.formStatus()));

        } else if (measurementDone.equalsIgnoreCase("yes")) {
            if (date == null) {
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "Response date and/or assessment date empty or in wrong format",
                        treatmentResponse.formStatus()));
            }
            if (response == null && !isFirstResponse) {
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "Measurement done is yes, but response is empty (non-first response)",
                        treatmentResponse.formStatus()));
            }
        } else if (measurementDone.equalsIgnoreCase("no")) {
            if (!treatmentResponse.isNotDoneResponse()) {
                if (date != null) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "Measurement done is no, but assessment date or response date is filled in",
                            treatmentResponse.formStatus(),
                            "Effective response date: " + date));
                }
                if (response != null) {
                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "Measurement done is no, but response filled in",
                            treatmentResponse.formStatus(),
                            "Response: " + response));
                }
            }
        } else {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Measurement done is not yes/no",
                    treatmentResponse.formStatus(),
                    "Measurement done: " + measurementDone));
        }

        if (response != null && date == null && !treatmentResponse.isNotDoneResponse()) {
            findings.add(ValidationFinding.of(ECRF_LEVEL,
                    patientIdentifier,
                    "Response filled in, but no assessment date and response date found",
                    treatmentResponse.formStatus(),
                    "Response: " + response));
        }
        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateDeathDate(@NotNull String patientIdentifier, @NotNull BaselineData baselineData,
            @NotNull List<BiopsyTreatmentData> treatments) {
        List<ValidationFinding> findings = Lists.newArrayList();
        LocalDate deathDate = baselineData.deathDate();
        if (deathDate != null && !treatments.isEmpty()) {
            treatments.sort(comparing(BiopsyTreatmentData::endDate, nullsLast(naturalOrder())));
            BiopsyTreatmentData lastTreatment = treatments.get(treatments.size() - 1);
            String treatmentGiven = lastTreatment.treatmentGiven();
            if (treatmentGiven != null && treatmentGiven.equalsIgnoreCase("yes")) {
                LocalDate lastTreatmentEndDate = lastTreatment.endDate();
                if (lastTreatmentEndDate == null || lastTreatmentEndDate.isAfter(deathDate)) {
                    String details = "Death date (" + deathDate + ") before end of last treatment (" + lastTreatmentEndDate + ")";

                    findings.add(ValidationFinding.of(ECRF_LEVEL,
                            patientIdentifier,
                            "Death date before end of last treatment",
                            FormStatus.merge(baselineData.deathStatus(), lastTreatment.formStatus()),
                            details));
                }
            }

        }
        return findings;
    }

    @NotNull
    @VisibleForTesting
    static List<ValidationFinding> validateInformedConsentDate(@NotNull String patientIdentifier, @NotNull BaselineData baselineData,
            @NotNull List<BiopsyData> biopsies) {
        List<ValidationFinding> findings = Lists.newArrayList();
        LocalDate informedConsentDate = baselineData.informedConsentDate();
        if (informedConsentDate != null && !biopsies.isEmpty()) {
            List<BiopsyData> biopsiesPriorToInformedConsent = biopsies.stream().filter(biopsy -> {
                LocalDate biopsyDate = biopsy.date();
                return biopsyDate != null && biopsyDate.isBefore(informedConsentDate);
            }).collect(Collectors.toList());
            if (biopsiesPriorToInformedConsent.size() > 0) {
                String detailsMessage =
                        "InformedConsentDate: " + informedConsentDate + ". biopsies: " + biopsiesPriorToInformedConsent.stream()
                                .map(BiopsyData::toString)
                                .collect(Collectors.toList());

                FormStatus mergedFormStatus = baselineData.informedConsentStatus();
                for (BiopsyData biopsy : biopsies) {
                    mergedFormStatus = FormStatus.merge(mergedFormStatus, biopsy.formStatus());
                }
                findings.add(ValidationFinding.of(ECRF_LEVEL,
                        patientIdentifier,
                        "At least 1 biopsy taken before informed consent date",
                        mergedFormStatus,
                        detailsMessage));
            }
        }
        return findings;
    }
}
