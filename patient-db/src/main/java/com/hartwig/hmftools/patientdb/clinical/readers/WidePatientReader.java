package com.hartwig.hmftools.patientdb.clinical.readers;

import java.time.LocalDate;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.consents.ConsentConfig;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedBiopsyType;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.clinical.datamodel.DrugData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableBiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableCuratedBiopsyType;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableDrugData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ImmutableValidationFinding;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.clinical.matchers.MatchResult;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideAvlTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideBiopsyData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideClinicalData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideEcrfModel;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideFiveDays;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WidePreAvlTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideResponseData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class WidePatientReader {

    private static final Logger LOGGER = LogManager.getLogger(WidePatientReader.class);

    @NotNull
    private final WideEcrfModel wideEcrfModel;
    @NotNull
    private final PrimaryTumorCurator primaryTumorCurator;
    @NotNull
    private final TreatmentCurator treatmentCurator;

    public WidePatientReader(@NotNull final WideEcrfModel wideEcrfModel, @NotNull final PrimaryTumorCurator primaryTumorCurator,
            @NotNull final TreatmentCurator treatmentCurator) {
        this.wideEcrfModel = wideEcrfModel;
        this.primaryTumorCurator = primaryTumorCurator;
        this.treatmentCurator = treatmentCurator;
    }

    @NotNull
    public Patient read(@NotNull String patientIdentifier, @Nullable String limsPrimaryTumorLocation,
            @NotNull List<SampleData> sequencedSamples, @NotNull Map<String, ConsentConfig> consentConfigMap) {
        BaselineData baseline = buildBaselineData(patientIdentifier, limsPrimaryTumorLocation, consentConfigMap);
        MatchResult<BiopsyData> matchedBiopsies = buildMatchedBiopsies(patientIdentifier, sequencedSamples);
        List<BiopsyTreatmentData> treatments = buildBiopsyTreatmentData(patientIdentifier);
        List<BiopsyTreatmentResponseData> responses = buildBiopsyTreatmentResponseData(patientIdentifier);

        return new Patient(patientIdentifier,
                baseline,
                ImmutablePreTreatmentData.builder().formStatus(FormStatus.undefined()).build(),
                sequencedSamples,
                matchedBiopsies.values(),
                treatments,
                responses,
                Lists.newArrayList(),
                Lists.newArrayList(),
                matchedBiopsies.findings());
    }

    @NotNull
    private BaselineData buildBaselineData(@NotNull String patientIdentifier, @Nullable String limsPrimaryTumorLocation,
            @NotNull Map<String, ConsentConfig> consentConfigMap) {
        List<WideFiveDays> fiveDays = dataForPatient(wideEcrfModel.fiveDays(), patientIdentifier);

        // There is one entry per biopsy but we assume they all share the same baseline information
        Integer birthYear = !fiveDays.isEmpty() ? fiveDays.get(0).birthYear() : null;
        String gender = !fiveDays.isEmpty() ? fiveDays.get(0).gender() : null;
        LocalDate informedConsentDate = !fiveDays.isEmpty() ? fiveDays.get(0).informedConsentDate() : null;

        CuratedPrimaryTumor curatedPrimaryTumor = primaryTumorCurator.search(patientIdentifier, limsPrimaryTumorLocation);

        ConsentConfig extractConsentConfigInfo = consentConfigMap.get("WIDE");

        LOGGER.info(extractConsentConfigInfo);
        return ImmutableBaselineData.builder()
                .registrationDate(null)
                .informedConsentDate(informedConsentDate)
                .pifVersion(extractConsentConfigInfo != null ? extractConsentConfigInfo.pifVersion() : null )
                .inDatabase(extractConsentConfigInfo != null ? extractConsentConfigInfo.inHMF() : null)
                .outsideEU(extractConsentConfigInfo != null ? extractConsentConfigInfo.outsideEU() : null)
                .gender(gender)
                .hospital(null)
                .birthYear(birthYear)
                .curatedPrimaryTumor(curatedPrimaryTumor)
                .deathDate(null)
                .demographyStatus(FormStatus.undefined())
                .primaryTumorStatus(FormStatus.undefined())
                .informedConsentStatus(FormStatus.undefined())
                .eligibilityStatus(FormStatus.undefined())
                .selectionCriteriaStatus(FormStatus.undefined())
                .deathStatus(FormStatus.undefined())
                .build();
    }

    @NotNull
    private MatchResult<BiopsyData> buildMatchedBiopsies(@NotNull String patientIdentifier, @NotNull List<SampleData> sequencedBiopsies) {
        List<WideBiopsyData> clinicalBiopsies = dataForPatient(wideEcrfModel.biopsies(), patientIdentifier);

        List<ValidationFinding> findings = Lists.newArrayList();
        if (clinicalBiopsies.size() < sequencedBiopsies.size()) {
            findings.add(biopsyMatchFinding(patientIdentifier,
                    "Not enough clinical biopsies to match for every sequenced sample",
                    "Clinical biopsies: " + clinicalBiopsies.size() + ", sequenced samples: " + sequencedBiopsies.size()));
        }

        List<BiopsyData> biopsies = Lists.newArrayList();
        List<WideFiveDays> fiveDaysForPatient = dataForPatient(wideEcrfModel.fiveDays(), patientIdentifier);
        for (WideBiopsyData clinicalBiopsy : clinicalBiopsies) {
            WideFiveDays fiveDays = lookupFiveDaysOnBiopsyDate(fiveDaysForPatient, clinicalBiopsy.biopsyDate());

            if (fiveDays != null) {
                BiopsyData finalBiopsy = toBiopsyData(fiveDays);
                SampleData sequencedSample = lookupOnPathologySampleId(sequencedBiopsies, clinicalBiopsy.pathologySampleId());
                if (sequencedSample != null && !sampleHasBeenMatched(biopsies, sequencedSample.sampleId())) {
                    finalBiopsy = ImmutableBiopsyData.builder().from(finalBiopsy).sampleId(sequencedSample.sampleId()).build();
                }

                biopsies.add(finalBiopsy);
            } else {
                findings.add(biopsyMatchFinding(patientIdentifier,
                        "Could not find five days entry for biopsy",
                        "WIDE Biopsy = " + clinicalBiopsy));
            }
        }

        return new MatchResult<>(biopsies, findings);
    }

    private static boolean sampleHasBeenMatched(@NotNull List<BiopsyData> biopsies, @NotNull String sampleId) {
        for (BiopsyData biopsy : biopsies) {
            String matchedSampleId = biopsy.sampleId();
            if (matchedSampleId != null && matchedSampleId.equals(sampleId)) {
                return true;
            }
        }
        return false;
    }

    @Nullable
    private static WideFiveDays lookupFiveDaysOnBiopsyDate(@NotNull List<WideFiveDays> fiveDaysForPatient, @Nullable LocalDate biopsyDate) {
        if (biopsyDate == null) {
            return null;
        }

        for (WideFiveDays fiveDays : fiveDaysForPatient) {
            LocalDate fiveDaysBiopsyDate = fiveDays.biopsyDate();
            if (fiveDaysBiopsyDate != null && fiveDaysBiopsyDate.equals(biopsyDate)) {
                return fiveDays;
            }
        }

        return null;
    }

    @Nullable
    private static SampleData lookupOnPathologySampleId(@NotNull List<SampleData> sequencedBiopsies, @NotNull String pathologySampleId) {
        String clinicalPathologySampleIdYear = extractYearFromPathologySampleId(pathologySampleId);
        String clinicalPathologySampleIdConvert = extractBiopsyIdFromPathologySampleId(pathologySampleId);

        if (clinicalPathologySampleIdYear != null && clinicalPathologySampleIdConvert != null) {
            for (SampleData sequencedBiopsy : sequencedBiopsies) {
                String limsPathologySampleIdYear = extractYearFromPathologySampleId(sequencedBiopsy.pathologySampleId());
                String limsPathologySampleIdConvert = extractBiopsyIdFromPathologySampleId(sequencedBiopsy.pathologySampleId());

                if (clinicalPathologySampleIdConvert.equals(limsPathologySampleIdConvert) && clinicalPathologySampleIdYear.equals(
                        limsPathologySampleIdYear)) {
                    return sequencedBiopsy;
                }
            }
        }
        return null;
    }

    @Nullable
    @VisibleForTesting
    static String extractYearFromPathologySampleId(@NotNull String pathologySampleId) {
        if (pathologySampleId.isEmpty()) {
            return null;
        }

        return pathologySampleId.split("-")[0];
    }

    @Nullable
    @VisibleForTesting
    static String extractBiopsyIdFromPathologySampleId(@NotNull String pathologySampleId) {
        if (pathologySampleId.isEmpty()) {
            return null;
        }

        String[] parts = pathologySampleId.split("-");
        if (parts.length >= 2) {
            return parts[1].replaceFirst("^0+(?!$)", "").split(" ")[0];
        } else {
            LOGGER.warn("Could not extract biopsy id from pathology sample id: {}", pathologySampleId);
            return null;
        }
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

    @NotNull
    private static BiopsyData toBiopsyData(@NotNull WideFiveDays fiveDays) {
        CuratedBiopsyType curatedBiopsyType = ImmutableCuratedBiopsyType.builder()
                .type("Unknown")
                .searchPrimaryTumorLocation(Strings.EMPTY)
                .searchCancerSubType(Strings.EMPTY)
                .searchBiopsySite(Strings.EMPTY)
                .searchBiopsyLocation(Strings.EMPTY)
                .build();

        return BiopsyData.of(fiveDays.biopsyDate(),
                null,
                null,
                curatedBiopsyType,
                fiveDays.biopsySite(),
                fiveDays.sampleTissue(),
                FormStatus.undefined());
    }

    @NotNull
    private List<BiopsyTreatmentData> buildBiopsyTreatmentData(@NotNull String patientIdentifier) {
        List<BiopsyTreatmentData> biopsyTreatmentDataList = Lists.newArrayList();

        for (WidePreAvlTreatmentData preAvlTreatment : dataForPatient(wideEcrfModel.preAvlTreatments(), patientIdentifier)) {
            biopsyTreatmentDataList.add(BiopsyTreatmentData.of(null,
                    "no",
                    null,
                    preAvlTreatmentDrugList(preAvlTreatment, treatmentCurator),
                    FormStatus.undefined()));
        }

        for (WideAvlTreatmentData avlTreatment : dataForPatient(wideEcrfModel.avlTreatments(), patientIdentifier)) {
            biopsyTreatmentDataList.add(BiopsyTreatmentData.of(null,
                    "yes",
                    null,
                    avlTreatmentDrugList(avlTreatment, treatmentCurator),
                    FormStatus.undefined()));
        }

        return biopsyTreatmentDataList;
    }

    @NotNull
    public static List<DrugData> preAvlTreatmentDrugList(@NotNull WidePreAvlTreatmentData preAvlTreatment,
            @NotNull TreatmentCurator treatmentCurator) {
        List<DrugData> drugs = Lists.newArrayList();

        LocalDate drugsEndDate = preAvlTreatment.lastSystemicTherapyDate();
        if (!preAvlTreatment.drug1().isEmpty()) {
            drugs.add(toSingleDrug(preAvlTreatment.drug1(), null, drugsEndDate, treatmentCurator));
        }

        if (!preAvlTreatment.drug2().isEmpty()) {
            drugs.add(toSingleDrug(preAvlTreatment.drug2(), null, drugsEndDate, treatmentCurator));
        }

        if (!preAvlTreatment.drug3().isEmpty()) {
            drugs.add(toSingleDrug(preAvlTreatment.drug3(), null, drugsEndDate, treatmentCurator));
        }

        if (!preAvlTreatment.drug4().isEmpty()) {
            drugs.add(toSingleDrug(preAvlTreatment.drug4(), null, drugsEndDate, treatmentCurator));
        }

        return drugs;
    }

    @NotNull
    private static List<DrugData> avlTreatmentDrugList(@NotNull WideAvlTreatmentData avlTreatment,
            @NotNull TreatmentCurator treatmentCurator) {
        List<DrugData> drugs = Lists.newArrayList();
        String drugName = avlTreatment.drug();

        if (!drugName.isEmpty()) {
            drugs.add(toSingleDrug(drugName, avlTreatment.startDate(), avlTreatment.endDate(), treatmentCurator));
        }
        return drugs;
    }

    @NotNull
    private static DrugData toSingleDrug(@NotNull String drug, @Nullable LocalDate startDate, @Nullable LocalDate endDate,
            @NotNull TreatmentCurator treatmentCurator) {
        return ImmutableDrugData.builder()
                .name(drug)
                .startDate(startDate)
                .endDate(endDate)
                .bestResponse(null)
                .curatedDrugs(treatmentCurator.search(drug))
                .build();
    }

    @NotNull
    private List<BiopsyTreatmentResponseData> buildBiopsyTreatmentResponseData(@NotNull String patientIdentifier) {
        List<BiopsyTreatmentResponseData> biopsyTreatmentResponseDataList = Lists.newArrayList();
        for (WideResponseData response : dataForPatient(wideEcrfModel.responses(), patientIdentifier)) {
            biopsyTreatmentResponseDataList.add(BiopsyTreatmentResponseData.of(null,
                    null,
                    response.date(),
                    determineResponse(response),
                    "yes",
                    null,
                    FormStatus.undefined()));
        }

        return biopsyTreatmentResponseDataList;
    }

    @NotNull
    @VisibleForTesting
    static String determineResponse(@NotNull WideResponseData response) {
        String responseString = Strings.EMPTY;
        if (response.recistDone()) {
            if (!response.recistResponse().isEmpty()) {
                responseString = "(" + response.timePoint() + ") " + response.recistResponse();
            }
        } else if (!response.noRecistResponse().isEmpty()) {
            responseString = "(" + response.timePoint() + ") " + response.noRecistResponse();
            if (!response.noRecistReasonStopTreatment().isEmpty()) {
                if (response.noRecistReasonStopTreatmentOther().isEmpty()) {
                    responseString = responseString + " (" + response.noRecistReasonStopTreatment() + ")";
                } else {
                    responseString = responseString + " (" + response.noRecistReasonStopTreatmentOther() + ")";
                }
            }
        }
        return responseString;
    }

    @NotNull
    private static <T extends WideClinicalData> List<T> dataForPatient(@NotNull List<T> data, @NotNull String patientIdentifier) {
        List<T> elementsForPatient = Lists.newArrayList();
        for (T element : data) {
            if (element.widePatientId().equals(patientIdentifier)) {
                elementsForPatient.add(element);
            }
        }
        return elementsForPatient;
    }
}