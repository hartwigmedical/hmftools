package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.CuratedBiopsyType;
import com.hartwig.hmftools.patientdb.data.CuratedDrug;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedBiopsyType;
import com.hartwig.hmftools.patientdb.data.ImmutableDrugData;
import com.hartwig.hmftools.patientdb.data.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.matchers.BiopsyMatcher;
import com.hartwig.hmftools.patientdb.matchers.MatchResult;
import com.hartwig.hmftools.patientdb.matchers.TreatmentMatcher;
import com.hartwig.hmftools.patientdb.matchers.TreatmentResponseMatcher;
import com.hartwig.hmftools.patientdb.readers.wide.WideAvlTreatmentData;
import com.hartwig.hmftools.patientdb.readers.wide.WideBiopsyData;
import com.hartwig.hmftools.patientdb.readers.wide.WideClinicalData;
import com.hartwig.hmftools.patientdb.readers.wide.WideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.WideFiveDays;
import com.hartwig.hmftools.patientdb.readers.wide.WidePreAvlTreatmentData;
import com.hartwig.hmftools.patientdb.readers.wide.WideResponseData;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class WidePatientReader {

    @NotNull
    private final WideEcrfModel wideEcrfModel;
    @NotNull
    private final TumorLocationCurator tumorLocationCurator;
    @NotNull
    private final TreatmentCurator treatmentCurator;

    public WidePatientReader(@NotNull final WideEcrfModel wideEcrfModel, @NotNull final TumorLocationCurator tumorLocationCurator,
            @NotNull final TreatmentCurator treatmentCurator) {
        this.wideEcrfModel = wideEcrfModel;
        this.tumorLocationCurator = tumorLocationCurator;
        this.treatmentCurator = treatmentCurator;
    }

    @NotNull
    public Patient read(@NotNull String patientIdentifier, @Nullable String limsPrimaryTumorLocation,
            @NotNull List<SampleData> sequencedSamples) {
        LocalDate biopsyDateCheck = null;
        for (WideBiopsyData biopsy : wideEcrfModel.biopsies()) {
            if (patientIdentifier.equals(biopsy.widePatientId())) {
                String limsPathologyTissueId = sequencedSamples.get(0).pathologySampleId();
                String limsPathologyTissueIdYear = extractYearFromPathologySampleId(limsPathologyTissueId);
                String limsPathologyTissueIdConvert = extractBiopsyIdFromPathologySampleId(limsPathologyTissueId);

                String ecrfPathologyTissueId = biopsy.pathologySampleId();
                String ecrfPathologyTissueIdYear = extractYearFromPathologySampleId(ecrfPathologyTissueId);
                String ecrfPathologyTissueIdConvert = extractBiopsyIdFromPathologySampleId(ecrfPathologyTissueId);

                if (ecrfPathologyTissueIdYear.equals(limsPathologyTissueIdYear) && ecrfPathologyTissueIdConvert.equals(
                        limsPathologyTissueIdConvert)) {
                    biopsyDateCheck = biopsy.biopsyDate();
                }
            }
        }

        LocalDate biopsyDate = bioptDate(wideEcrfModel.biopsies(), patientIdentifier);

        List<BiopsyData> biopsyData = toBiopsyData(wideEcrfModel.fiveDays(), patientIdentifier, biopsyDateCheck);

        MatchResult<BiopsyData> matchedBiopsies =
                BiopsyMatcher.matchBiopsiesToTumorSamples(patientIdentifier, sequencedSamples, biopsyData);

        MatchResult<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies(patientIdentifier,
                withSampleMatchOnly(matchedBiopsies),
                toBiopsyTreatmentData(wideEcrfModel.avlTreatments(), treatmentCurator, patientIdentifier, biopsyDate));

        // We also match responses to unmatched treatments. Not sure that is optimal. See also DEV-477.
        MatchResult<BiopsyTreatmentResponseData> matchedResponses = TreatmentResponseMatcher.matchTreatmentResponsesToTreatments(
                patientIdentifier,
                matchedTreatments.values(),
                toBiopsyTreatmentResponseData(wideEcrfModel.responses(), patientIdentifier));

        List<ValidationFinding> findings = Lists.newArrayList();
        findings.addAll(matchedBiopsies.findings());
        findings.addAll(matchedTreatments.findings());
        findings.addAll(matchedResponses.findings());

        return new Patient(patientIdentifier,
                buildBaselineData(patientIdentifier, limsPrimaryTumorLocation),
                preTreatmentData(wideEcrfModel.preAvlTreatments(),
                        treatmentCurator,
                        patientIdentifier,
                        biopsyDate,
                        wideEcrfModel.avlTreatments()),
                sequencedSamples,
                matchedBiopsies.values(),
                matchedTreatments.values(),
                matchedResponses.values(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                findings);
    }

    @NotNull
    @VisibleForTesting
    static String extractYearFromPathologySampleId(@NotNull String pathologySampleId) {
        return pathologySampleId.split("-")[0];
    }

    @NotNull
    @VisibleForTesting
    static String extractBiopsyIdFromPathologySampleId(@NotNull String pathologySampleId) {
        return pathologySampleId.split("-")[1].replaceFirst("^0+(?!$)", "").split(" ")[0];
    }

    @NotNull
    private static List<BiopsyData> withSampleMatchOnly(@NotNull MatchResult<BiopsyData> biopsies) {
        List<BiopsyData> biopsiesWithMatchedSample = Lists.newArrayList();
        for (BiopsyData biopsy : biopsies.values()) {
            if (biopsy.sampleId() != null) {
                biopsiesWithMatchedSample.add(biopsy);
            }
        }
        return biopsiesWithMatchedSample;
    }

    @NotNull
    private BaselineData buildBaselineData(@NotNull String patientIdentifier, @Nullable String limsPrimaryTumorLocation) {
        CuratedTumorLocation curatedTumorLocation = tumorLocationCurator.search(limsPrimaryTumorLocation);

        Integer birthYear = null;
        String gender = Strings.EMPTY;
        LocalDate informedConsentDate = null;

        for (WideFiveDays fiveDays : wideEcrfModel.fiveDays()) {
            if (patientIdentifier.equals(fiveDays.widePatientId())) {
                birthYear = fiveDays.birthYear();
                gender = fiveDays.gender();
                informedConsentDate = fiveDays.informedConsentDate();
            }
        }

        return ImmutableBaselineData.builder()
                .registrationDate(null)
                .informedConsentDate(informedConsentDate)
                .gender(gender)
                .hospital(null)
                .birthYear(birthYear)
                .curatedTumorLocation(curatedTumorLocation)
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
    private MatchResult<BiopsyData> matchedBiopsies(@NotNull String patientIdentifier, @NotNull List<SampleData> sequencedSamples) {
        List<WideFiveDays> fiveDaysForPatient = dataForPatient(wideEcrfModel.fiveDays(), patientIdentifier);
        List<WideBiopsyData> biopsiesForPatient = dataForPatient(wideEcrfModel.biopsies(), patientIdentifier);

        List<ValidationFinding> findings = Lists.newArrayList();

        LocalDate biopsyDateCheck = null;
        for (WideBiopsyData biopsy : wideEcrfModel.biopsies()) {
            if (patientIdentifier.equals(biopsy.widePatientId())) {
                String limsPathologyTissueId = sequencedSamples.get(0).pathologySampleId();
                String limsPathologyTissueIdYear = extractYearFromPathologySampleId(limsPathologyTissueId);
                String limsPathologyTissueIdConvert = extractBiopsyIdFromPathologySampleId(limsPathologyTissueId);

                String ecrfPathologyTissueId = biopsy.pathologySampleId();
                String ecrfPathologyTissueIdYear = extractYearFromPathologySampleId(ecrfPathologyTissueId);
                String ecrfPathologyTissueIdConvert = extractBiopsyIdFromPathologySampleId(ecrfPathologyTissueId);

                if (ecrfPathologyTissueIdYear.equals(limsPathologyTissueIdYear) && ecrfPathologyTissueIdConvert.equals(
                        limsPathologyTissueIdConvert)) {
                    biopsyDateCheck = biopsy.biopsyDate();
                }
            }
        }

        LocalDate biopsyDate = bioptDate(wideEcrfModel.biopsies(), patientIdentifier);

        List<BiopsyData> biopsyData = toBiopsyData(wideEcrfModel.fiveDays(), patientIdentifier, biopsyDateCheck);

        MatchResult<BiopsyData> matchedBiopsies =
                BiopsyMatcher.matchBiopsiesToTumorSamples(patientIdentifier, sequencedSamples, biopsyData);

        return null;
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

    @NotNull
    private static PreTreatmentData preTreatmentData(@NotNull List<WidePreAvlTreatmentData> widePreAvlTreatmentDatum2s,
            @NotNull TreatmentCurator treatmentCurator, @NotNull String patientIdentifier, @NotNull LocalDate biopsyDate,
            @NotNull List<WideAvlTreatmentData> treatmentData) {
        return ImmutablePreTreatmentData.of(null,
                null,
                readDrugsPreTreatment(widePreAvlTreatmentDatum2s, treatmentCurator, patientIdentifier, biopsyDate, treatmentData),
                FormStatus.undefined());
    }

    @NotNull
    public static List<DrugData> readDrugsPreTreatment(@NotNull List<WidePreAvlTreatmentData> preTreatmentData,
            @NotNull TreatmentCurator treatmentCurator, @NotNull String patientIdentifier, @NotNull LocalDate biopsyDate,
            @NotNull List<WideAvlTreatmentData> treatmentData) {
        List<DrugData> drugs = Lists.newArrayList();
        for (WidePreAvlTreatmentData preTreatment : preTreatmentData) {
            if (patientIdentifier.equals(preTreatment.widePatientId())) {
                LocalDate drugsEndDate = preTreatment.lastSystemicTherapyDate();
                String drugsName = Strings.EMPTY;

                if (!preTreatment.drug1().equals(Strings.EMPTY)) {
                    if (!drugsName.isEmpty()) {
                        drugsName = drugsName + "," + preTreatment.drug1();
                    } else {
                        drugsName = preTreatment.drug1();
                    }
                }

                if (!preTreatment.drug2().equals(Strings.EMPTY)) {
                    if (!drugsName.isEmpty()) {
                        drugsName = drugsName + "," + preTreatment.drug2();
                    } else {
                        drugsName = preTreatment.drug2();
                    }
                }

                if (!preTreatment.drug3().equals(Strings.EMPTY)) {
                    if (!drugsName.isEmpty()) {
                        drugsName = drugsName + "," + preTreatment.drug3();
                    } else {
                        drugsName = preTreatment.drug3();
                    }
                }

                if (!preTreatment.drug4().equals(Strings.EMPTY)) {
                    if (!drugsName.isEmpty()) {
                        drugsName = drugsName + "," + preTreatment.drug4();
                    } else {
                        drugsName = preTreatment.drug4();
                    }

                }
                List<CuratedDrug> curatedDrugs = treatmentCurator.search(drugsName);
                drugs.add(ImmutableDrugData.of(drugsName, null, drugsEndDate, null, curatedDrugs));
            }
        }

        for (WideAvlTreatmentData postTreatment : treatmentData) {
            if (patientIdentifier.equals(postTreatment.widePatientId())) {
                if (postTreatment.startDate().isBefore(biopsyDate)) {
                    List<CuratedDrug> curatedDrugs = treatmentCurator.search(postTreatment.drug());
                    drugs.add(ImmutableDrugData.of(postTreatment.drug(),
                            postTreatment.startDate(),
                            postTreatment.endDate(),
                            null,
                            curatedDrugs));
                }
            }
        }
        return drugs;
    }

    @NotNull
    private static LocalDate bioptDate(@NotNull List<WideBiopsyData> wideBiopsyData, @NotNull String patientIdentifier) {
        LocalDate biopsyDate = LocalDate.now();
        for (WideBiopsyData biopsy : wideBiopsyData) {
            if (patientIdentifier.equals(biopsy.widePatientId())) {
                biopsyDate = biopsy.biopsyDate();
            }
        }
        return biopsyDate;
    }

    @NotNull
    private static List<BiopsyData> toBiopsyData(@NotNull List<WideFiveDays> wideBiopsyData, @NotNull String patientIdentifier,
            @Nullable LocalDate biopsyCheckDate) {
        List<BiopsyData> biopsyDataList = Lists.newArrayList();
        CuratedBiopsyType curatedBiopsyType = ImmutableCuratedBiopsyType.builder()
                .type("Unknown")
                .searchPrimaryTumorLocation(Strings.EMPTY)
                .searchCancerSubType(Strings.EMPTY)
                .searchBiopsySite(Strings.EMPTY)
                .searchBiopsyLocation(Strings.EMPTY)
                .build();

        for (WideFiveDays biopsyData : wideBiopsyData) {
            if (patientIdentifier.equals(biopsyData.widePatientId())) {
                if (biopsyData.biopsyDate() != null && biopsyCheckDate != null) {
                    if (biopsyData.biopsyDate().equals(biopsyCheckDate)) {
                        biopsyDataList.add(BiopsyData.of(biopsyData.biopsyDate(),
                                null,
                                null,
                                curatedBiopsyType,
                                biopsyData.biopsySite(),
                                biopsyData.sampleTissue(),
                                FormStatus.undefined()));
                    }
                }
            }
        }
        return biopsyDataList;

    }

    @NotNull
    private static List<BiopsyTreatmentData> toBiopsyTreatmentData(@NotNull List<WideAvlTreatmentData> wideAvlTreatmentData,
            @NotNull TreatmentCurator treatmentCurator, @NotNull String patientIdentifier, @NotNull LocalDate biopsyDate) {
        List<BiopsyTreatmentData> biopsyTreatmentDataList = Lists.newArrayList();
        for (WideAvlTreatmentData treatmentData : wideAvlTreatmentData) {
            if (patientIdentifier.equals(treatmentData.widePatientId())) {
                if (biopsyDate.isEqual((treatmentData.startDate())) || treatmentData.startDate().isAfter(biopsyDate)) {
                    biopsyTreatmentDataList.add(BiopsyTreatmentData.of(null,
                            "yes",
                            null,
                            readDrugsPostTreatment(treatmentData, treatmentCurator),
                            FormStatus.undefined()));
                }
            }
        }
        return biopsyTreatmentDataList;
    }

    @NotNull
    private static List<DrugData> readDrugsPostTreatment(@NotNull WideAvlTreatmentData treatmentData,
            @NotNull TreatmentCurator treatmentCurator) {
        List<DrugData> drugs = Lists.newArrayList();
        String drugName = treatmentData.drug();
        LocalDate drugsStartDate = treatmentData.startDate();
        LocalDate drugsEndDate = treatmentData.endDate();

        if (!drugName.isEmpty() || drugsStartDate != null || drugsEndDate != null) {
            List<CuratedDrug> curatedDrugs = drugName.isEmpty() ? Lists.newArrayList() : treatmentCurator.search(drugName);
            drugs.add(ImmutableDrugData.of(drugName, drugsStartDate, drugsEndDate, null, curatedDrugs));
        }
        return drugs;
    }

    @NotNull
    private static List<BiopsyTreatmentResponseData> toBiopsyTreatmentResponseData(@NotNull List<WideResponseData> wideResponseData,
            @NotNull String patientIdentifier) {
        List<BiopsyTreatmentResponseData> biopsyTreatmentResponseDataList = Lists.newArrayList();
        for (WideResponseData responseData : wideResponseData) {
            if (patientIdentifier.equals(responseData.widePatientId())) {
                biopsyTreatmentResponseDataList.add(BiopsyTreatmentResponseData.of(null,
                        null,
                        responseData.date(),
                        determineResponse(responseData),
                        "yes",
                        null,
                        FormStatus.undefined()));
            }
        }

        return biopsyTreatmentResponseDataList;
    }

    @NotNull
    public static String determineResponse(@NotNull WideResponseData responseData) {
        String response;
        if (responseData.recistDone()) {
            response = "(" + responseData.timePoint() + ") " + responseData.recistResponse();
        } else {
            response = "(" + responseData.timePoint() + ") " + responseData.noRecistResponse();
            if (!response.contains("clinical benefit (continuation of treatment)") && !response.isEmpty()) {
                response = response + " (" + responseData.noRecistReasonStopTreatment() + ")";
                if (!responseData.noRecistReasonStopTreatmentOther().isEmpty()) {
                    response = response.replace(")", " + " + responseData.noRecistReasonStopTreatmentOther() + ")");
                }
            }
        }
        return response;
    }
}