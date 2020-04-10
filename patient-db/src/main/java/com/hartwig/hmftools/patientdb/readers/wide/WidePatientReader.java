package com.hartwig.hmftools.patientdb.readers.wide;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeFormatterBuilder;
import java.util.List;
import java.util.Locale;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.curators.BiopsySiteCurator;
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
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.ImmutableDrugData;
import com.hartwig.hmftools.patientdb.data.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.matchers.BiopsyMatcher;
import com.hartwig.hmftools.patientdb.matchers.MatchResult;
import com.hartwig.hmftools.patientdb.matchers.TreatmentMatcher;
import com.hartwig.hmftools.patientdb.matchers.TreatmentResponseMatcher;

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
    private final TumorLocationCurator tumorLocationCurator;
    @NotNull
    private final BiopsySiteCurator biopsySiteCurator;
    @NotNull
    private final TreatmentCurator treatmentCurator;

    public WidePatientReader(@NotNull final WideEcrfModel wideEcrfModel, @NotNull final TumorLocationCurator tumorLocationCurator,
            @NotNull BiopsySiteCurator biopsySiteCurator, @NotNull final TreatmentCurator treatmentCurator) {
        this.wideEcrfModel = wideEcrfModel;
        this.tumorLocationCurator = tumorLocationCurator;
        this.biopsySiteCurator = biopsySiteCurator;
        this.treatmentCurator = treatmentCurator;
    }

    @NotNull
    public static String extractYearOfTissueId(@NotNull String tissueId) {
        return tissueId.split("-")[0];
    }

    @NotNull
    public static String extractBiopsyIdOfTissueId(@NotNull String tissueId) {
        return tissueId.split("-")[1].replaceFirst("^0+(?!$)", "").split(" ")[0];

    }

    @NotNull
    public Patient read(@NotNull String patientIdentifier, @Nullable String primaryTumorLocation,
            @NotNull List<SampleData> sequencedSamples, @NotNull Lims lims, @NotNull String tumorBarcode) {

        String biopsyDateCheck = Strings.EMPTY;
        for (WideBiopsyData biopsy : wideEcrfModel.biopsies()) {
            if (patientIdentifier.equals(biopsy.wideId())) {
                String limsPathologyTissueId = lims.hospitalPathologySampleId(tumorBarcode);
                String limsPathologyTissueIdYear = extractYearOfTissueId(limsPathologyTissueId);
                String limsPathologyTissueIdConvert = extractBiopsyIdOfTissueId(limsPathologyTissueId);

                String ecrfPathologyTissueId = biopsy.tissueId();
                String ecrfPathologyTissueIdYear = extractYearOfTissueId(ecrfPathologyTissueId);
                String ecrfPathologyTissueIdConvert = extractBiopsyIdOfTissueId(ecrfPathologyTissueId);

                if (ecrfPathologyTissueIdYear.equals(limsPathologyTissueIdYear) && ecrfPathologyTissueIdConvert.equals(
                        limsPathologyTissueIdConvert)) {
                    biopsyDateCheck = createInterpretDateNL(biopsy.bioptDate()).toString();
                }
            }
        }

        String birthYear = Strings.EMPTY;
        String gender = Strings.EMPTY;
        String informedConsent = Strings.EMPTY;
        String sharedData = Strings.EMPTY;
        String biopsySite = Strings.EMPTY;
        String sampleTissue = Strings.EMPTY;

        for (WideFiveDays fiveDays : wideEcrfModel.fiveDays()) {
            if (patientIdentifier.equals(fiveDays.wideId())) {
                if (!fiveDays.dateBiopsy().equals(Strings.EMPTY) && !biopsyDateCheck.equals(Strings.EMPTY)) {
                    if (createInterpretDateEN(fiveDays.dateBiopsy()).equals(convertStringToLocalDate(biopsyDateCheck))) {
                        birthYear = fiveDays.birthYear();
                        gender = fiveDays.gender();
                        informedConsent = fiveDays.informedConsent();
                        sharedData = fiveDays.dateShared();
                        biopsySite = fiveDays.biopsySite();
                        sampleTissue = fiveDays.sampleTissue();
                    }
                }
            }
        }

        LocalDate biopsyDate = bioptDate(wideEcrfModel.biopsies(), patientIdentifier);

        List<BiopsyData> biopsyData = toBiopsyData(wideEcrfModel.fiveDays(),
                biopsySiteCurator,
                patientIdentifier,
                biopsyDateCheck,
                biopsySite,
                sampleTissue,
                tumorLocationCurator.search(primaryTumorLocation));

        MatchResult<BiopsyData> matchedBiopsies =
                BiopsyMatcher.matchBiopsiesToTumorSamples(patientIdentifier, sequencedSamples, biopsyData);

        MatchResult<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies(patientIdentifier,
                withSampleMatchOnly(matchedBiopsies),
                toBiopsyTreatmentData(wideEcrfModel.treatments(), treatmentCurator, patientIdentifier, biopsyDate));

        // We also match responses to unmatched treatments. Not sure that is optimal. See also DEV-477.
        MatchResult<BiopsyTreatmentResponseData> matchedResponses = TreatmentResponseMatcher.matchTreatmentResponsesToTreatments(
                patientIdentifier,
                matchedTreatments.values(),
                toBiopsyTreatmentResponseData(wideEcrfModel.responses(), patientIdentifier));

        final List<ValidationFinding> findings = Lists.newArrayList();
        //        findings.addAll(matchedBiopsies.findings());
        //        findings.addAll(matchedTreatments.findings());
        //        findings.addAll(matchedResponses.findings());

        return new Patient(patientIdentifier,
                toBaselineData(tumorLocationCurator.search(primaryTumorLocation), birthYear, gender, informedConsent, sharedData),
                preTreatmentData(wideEcrfModel.preTreatments(),
                        treatmentCurator,
                        patientIdentifier,
                        biopsyDate,
                        wideEcrfModel.treatments()),
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
    public static LocalDate createInterpretDateIC(@NotNull String date) {
        String format = "dd/MM/yyyy";
        return LocalDate.parse(date, DateTimeFormatter.ofPattern(format));
    }

    @NotNull
    @VisibleForTesting
    public static LocalDate convertStringToLocalDate(@NotNull String date) {
        String format = "yyyy-MM-dd";
        return LocalDate.parse(date, DateTimeFormatter.ofPattern(format));
    }

    @NotNull
    @VisibleForTesting
    public static LocalDate createInterpretDateNL(@NotNull String date) {
        DateTimeFormatter inputFormatter =
                new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(new Locale("nl", "NL"));
        LocalDate localDate = LocalDate.parse(date, inputFormatter);

        DateTimeFormatter outputFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
        String formattedString = localDate.format(outputFormatter);
        return LocalDate.parse(formattedString);
    }

    @NotNull
    @VisibleForTesting
    public static LocalDate createInterpretDateEN(@NotNull String date) {
        DateTimeFormatter inputFormatter =
                new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(Locale.ENGLISH);
        LocalDate localDate = LocalDate.parse(date, inputFormatter);

        DateTimeFormatter outputFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
        String formattedString = localDate.format(outputFormatter);
        return LocalDate.parse(formattedString);
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
    private static BaselineData toBaselineData(@NotNull CuratedTumorLocation curatedTumorLocation, @NotNull String birthYear,
            @NotNull String gender, @NotNull String informedConsent, @NotNull String sharedData) {
        return ImmutableBaselineData.of(null,
                sharedData.equals("1") && !informedConsent.isEmpty() ? createInterpretDateIC(informedConsent) : null,
                convertGender(gender),
                null,
                birthYear.isEmpty() ? null : Integer.valueOf(birthYear),
                curatedTumorLocation,
                null,
                FormStatus.undefined(),
                FormStatus.undefined(),
                FormStatus.undefined(),
                FormStatus.undefined(),
                FormStatus.undefined(),
                FormStatus.undefined());
    }

    @NotNull
    private static String convertGender(@NotNull String gender) {
        if (gender.equals("1")) {
            return "Male";
        } else if (gender.equals("2")) {
            return "Female";
        } else {
            return Strings.EMPTY;
        }
    }

    @NotNull
    private static PreTreatmentData preTreatmentData(@NotNull List<WidePreTreatmentData> widePreTreatmentData,
            @NotNull final TreatmentCurator treatmentCurator, @NotNull String patientIdentifier, @NotNull LocalDate biopsyDate,
            @NotNull List<WideTreatmentData> treatmentData) {

        return ImmutablePreTreatmentData.of(null,
                null,
                readDrugsPreTreatment(widePreTreatmentData, treatmentCurator, patientIdentifier, biopsyDate, treatmentData),
                FormStatus.undefined());
    }

    @NotNull
    private static List<DrugData> readDrugsPreTreatment(@NotNull List<WidePreTreatmentData> preTreatmentData,
            @NotNull final TreatmentCurator treatmentCurator, @NotNull String patientIdentifier, @NotNull LocalDate biopsyDate,
            @NotNull List<WideTreatmentData> treatmentData) {
        List<DrugData> drugs = Lists.newArrayList();
        for (WidePreTreatmentData preTreatment : preTreatmentData) {
            if (patientIdentifier.equals(preTreatment.patientId())) {
                LocalDate drugsEndDate = preTreatment.dateLastSystemicTherapy().isEmpty()
                        ? null
                        : createInterpretDateNL(preTreatment.dateLastSystemicTherapy());
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
                final List<CuratedDrug> curatedDrugs = treatmentCurator.search(drugsName);
                drugs.add(ImmutableDrugData.of(drugsName, null, drugsEndDate, null, curatedDrugs));

            }
        }

        for (WideTreatmentData postTreatment : treatmentData) {
            if (patientIdentifier.equals(postTreatment.sampleId())) {
                if (createInterpretDateNL(postTreatment.startDate()).isBefore(biopsyDate)) {
                    final List<CuratedDrug> curatedDrugs = treatmentCurator.search(postTreatment.drug());
                    drugs.add(ImmutableDrugData.of(postTreatment.drug(),
                            postTreatment.startDate().equals(Strings.EMPTY) ? null : createInterpretDateNL(postTreatment.startDate()),
                            postTreatment.endDate().equals(Strings.EMPTY) ? null : createInterpretDateNL(postTreatment.endDate()),
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
            if (patientIdentifier.equals(biopsy.wideId())) {
                biopsyDate = createInterpretDateNL(biopsy.bioptDate());
            }
        }
        return biopsyDate;
    }

    @NotNull
    private static List<BiopsyData> toBiopsyData(@NotNull List<WideFiveDays> wideBiopsyData, @NotNull BiopsySiteCurator biopsySiteCurator,
            @NotNull String patientIdentifier, @NotNull String biopsyCheck, @NotNull String biopsySite, @NotNull String sampleTissue,
            @NotNull CuratedTumorLocation curatedTumorLocation) {
        List<BiopsyData> biopsyDataList = Lists.newArrayList();
        final CuratedBiopsyType curatedBiopsyType = biopsySiteCurator.search(curatedTumorLocation.primaryTumorLocation(),
                curatedTumorLocation.subType(),
                biopsySite,
                sampleTissue);

        for (WideFiveDays biopsyData : wideBiopsyData) {
            if (patientIdentifier.equals(biopsyData.wideId())) {
                if (!biopsyData.dateBiopsy().equals(Strings.EMPTY) && !biopsyCheck.equals(Strings.EMPTY)) {
                    if (createInterpretDateEN(biopsyData.dateBiopsy()).equals(convertStringToLocalDate(biopsyCheck))) {
                        biopsyDataList.add(ImmutableBiopsyData.of(biopsyData.dateBiopsy().isEmpty()
                                        ? null
                                        : createInterpretDateEN(biopsyData.dateBiopsy()),
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
        LOGGER.info(biopsyDataList);
        return biopsyDataList;

    }

    @NotNull
    private static List<BiopsyTreatmentData> toBiopsyTreatmentData(@NotNull List<WideTreatmentData> wideTreatmentData,
            @NotNull final TreatmentCurator treatmentCurator, @NotNull String patientIdentifier, @NotNull LocalDate biopsyDate) {
        List<BiopsyTreatmentData> biopsyTreatmentDataList = Lists.newArrayList();
        for (WideTreatmentData treatmentData : wideTreatmentData) {
            if (patientIdentifier.equals(treatmentData.sampleId())) {
                if (biopsyDate.isEqual(createInterpretDateNL(treatmentData.startDate())) || createInterpretDateNL(treatmentData.startDate())
                        .isAfter(biopsyDate)) {
                    biopsyTreatmentDataList.add(BiopsyTreatmentData.of("yes",
                            null,
                            readDrugsPostTreatment(treatmentData, treatmentCurator),
                            FormStatus.undefined()));
                }
            }
        }
        LOGGER.info(biopsyTreatmentDataList);
        return biopsyTreatmentDataList;
    }

    @NotNull
    private static List<DrugData> readDrugsPostTreatment(@NotNull WideTreatmentData treatmentData,
            @NotNull final TreatmentCurator treatmentCurator) {
        final List<DrugData> drugs = Lists.newArrayList();
        String drugName = treatmentData.drug();
        LocalDate drugsStartDate = treatmentData.startDate().isEmpty() ? null : createInterpretDateNL(treatmentData.startDate());
        LocalDate drugsEndDate = treatmentData.endDate().isEmpty() ? null : createInterpretDateNL(treatmentData.endDate());

        if (!drugName.isEmpty() || drugsStartDate != null || drugsEndDate != null) {
            final List<CuratedDrug> curatedDrugs = drugName.isEmpty() ? Lists.newArrayList() : treatmentCurator.search(drugName);
            drugs.add(ImmutableDrugData.of(drugName, drugsStartDate, drugsEndDate, null, curatedDrugs));
        }
        return drugs;
    }

    @NotNull
    private static List<BiopsyTreatmentResponseData> toBiopsyTreatmentResponseData(@NotNull List<WideResponseData> wideResponseData,
            @NotNull String patientIdentifier) {
        List<BiopsyTreatmentResponseData> biopsyTreatmentResponseDataList = Lists.newArrayList();
        for (WideResponseData responseData : wideResponseData) {
            if (patientIdentifier.equals(responseData.patientId())) {
                biopsyTreatmentResponseDataList.add(ImmutableBiopsyTreatmentResponseData.of(null,
                        createInterpretDateNL(responseData.date().isEmpty() ? Strings.EMPTY : responseData.date()),
                        determineResponse(responseData),
                        "yes",
                        null,
                        FormStatus.undefined()));
            }
        }

        LOGGER.info(biopsyTreatmentResponseDataList);
        return biopsyTreatmentResponseDataList;
    }

    @NotNull
    private static String determineResponse(@NotNull WideResponseData responseData) {
        String response = Strings.EMPTY;
        if (responseData.recistNotDone().equals("FALSE")) {
            response = responseData.responseAccordingRecist();
        } else if (responseData.recistNotDone().equals("WAAR")) { //change to TRUE
            response = responseData.clinicalDecision();
            if (!response.contains("clinical benifit (continuation of treatment)") && !response.isEmpty()) {
                response = response + " (" + responseData.reasonStopTreatment() + ")";
                if (!responseData.reasonStopTreatmentOther().isEmpty()) {
                    response = response.replace(")", " + " + responseData.reasonStopTreatmentOther() + ")");
                }
            }
        }
        return response;
    }
}