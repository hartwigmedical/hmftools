package com.hartwig.hmftools.patientdb.readers.wide;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeFormatterBuilder;
import java.util.List;
import java.util.Locale;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
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
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;
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

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class WidePatientReader {

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
    public Patient read(@NotNull String patientIdentifier, @Nullable String primaryTumorLocation,
            @NotNull List<SampleData> sequencedSamples) {
        final MatchResult<BiopsyData> matchedBiopsies = BiopsyMatcher.matchBiopsiesToTumorSamples(patientIdentifier,
                sequencedSamples,
                toBiopsyData(wideEcrfModel.biopsies(), biopsySiteCurator, patientIdentifier));

        final MatchResult<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies(patientIdentifier,
                withSampleMatchOnly(matchedBiopsies),
                toBiopsyTreatmentData(wideEcrfModel.treatments(), treatmentCurator, patientIdentifier));

        // We also match responses to unmatched treatments. Not sure that is optimal. See also DEV-477.
        final MatchResult<BiopsyTreatmentResponseData> matchedResponses = TreatmentResponseMatcher.matchTreatmentResponsesToTreatments(
                patientIdentifier,
                matchedTreatments.values(),
                toBiopsyTreatmentResponseData(wideEcrfModel.responses(), patientIdentifier));

        return new Patient(patientIdentifier,
                toBaselineData(tumorLocationCurator.search(primaryTumorLocation)),
                preTreatmentData(wideEcrfModel.preTreatments(), treatmentCurator, patientIdentifier),
                sequencedSamples,
                matchedBiopsies.values(),
                matchedTreatments.values(),
                matchedResponses.values(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList());
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
    private static BaselineData toBaselineData(@NotNull CuratedTumorLocation curatedTumorLocation) {
        return ImmutableBaselineData.of(null,
                null,
                null,
                null,
                null,
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
    private static PreTreatmentData preTreatmentData(@NotNull List<WidePreTreatmentData> widePreTreatmentData,
            @NotNull final TreatmentCurator treatmentCurator, @NotNull String patientIdentifier) {
        return ImmutablePreTreatmentData.of(null,
                null,
                readDrugsPreTreatment(widePreTreatmentData, treatmentCurator, patientIdentifier),
                FormStatus.undefined());
    }

    @NotNull
    private static List<DrugData> readDrugsPreTreatment(@NotNull List<WidePreTreatmentData> preTreatmentData,
            @NotNull final TreatmentCurator treatmentCurator, @NotNull String patientIdentifier) {
        List<DrugData> drugs = Lists.newArrayList();
        for (WidePreTreatmentData preTreatment : preTreatmentData) {
            if (patientIdentifier.equals(preTreatment.patientId())) {
                LocalDate drugsEndDate = preTreatment.dateLastSystemicTherapy().isEmpty()
                        ? null
                        : createInterpretDate(preTreatment.dateLastSystemicTherapy());

                if (!preTreatment.drug1().isEmpty() || drugsEndDate != null) {
                    final List<CuratedDrug> curatedDrugs1 = treatmentCurator.search(preTreatment.drug1());
                    drugs.add(ImmutableDrugData.of(preTreatment.drug1(), null, drugsEndDate, null, curatedDrugs1));
                }

                if (!preTreatment.drug2().isEmpty() || drugsEndDate != null) {
                    final List<CuratedDrug> curatedDrugs2 = treatmentCurator.search(preTreatment.drug2());
                    drugs.add(ImmutableDrugData.of(preTreatment.drug2(), null, drugsEndDate, null, curatedDrugs2));
                }

                if (!preTreatment.drug3().isEmpty() || drugsEndDate != null) {
                    final List<CuratedDrug> curatedDrugs3 = treatmentCurator.search(preTreatment.drug3());
                    drugs.add(ImmutableDrugData.of(preTreatment.drug3(), null, drugsEndDate, null, curatedDrugs3));
                }

                if (!preTreatment.drug4().isEmpty() || drugsEndDate != null) {
                    final List<CuratedDrug> curatedDrugs4 = treatmentCurator.search(preTreatment.drug4());
                    drugs.add(ImmutableDrugData.of(preTreatment.drug4(), null, drugsEndDate, null, curatedDrugs4));
                }
            }
        }
        return drugs;
    }

    @NotNull
    @VisibleForTesting
    public static LocalDate createInterpretDate(@NotNull String date) {
        DateTimeFormatter inputFormatter =
                new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(new Locale("nl", "NL"));
        LocalDate localDate = LocalDate.parse(date, inputFormatter);

        DateTimeFormatter outputFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
        String formattedString = localDate.format(outputFormatter);
        return LocalDate.parse(formattedString);
    }

    @NotNull
    private static List<BiopsyData> toBiopsyData(@NotNull List<WideBiopsyData> wideBiopsyData, @NotNull BiopsySiteCurator biopsySiteCurator,
            @NotNull String patientIdentifier) {
        List<BiopsyData> biopsyDataList = Lists.newArrayList();
        final CuratedBiopsyType curatedBiopsyType = biopsySiteCurator.search(Strings.EMPTY, Strings.EMPTY, Strings.EMPTY, Strings.EMPTY);

        for (WideBiopsyData biopsyData : wideBiopsyData) {
            if (patientIdentifier.equals(biopsyData.wideId())) {
                biopsyDataList.add(ImmutableBiopsyData.of(biopsyData.bioptDate().isEmpty()
                        ? null
                        : createInterpretDate(biopsyData.bioptDate()), null, null, curatedBiopsyType, null, null, FormStatus.undefined()));
            }
        }
        return biopsyDataList;

    }

    @NotNull
    private static List<BiopsyTreatmentData> toBiopsyTreatmentData(@NotNull List<WideTreatmentData> wideTreatmentData,
            @NotNull final TreatmentCurator treatmentCurator, @NotNull String patientIdentifier) {
        List<BiopsyTreatmentData> biopsyTreatmentDataList = Lists.newArrayList();
        for (WideTreatmentData treatmentData : wideTreatmentData) {
            if (patientIdentifier.equals(treatmentData.sampleId())) {
                biopsyTreatmentDataList.add(ImmutableBiopsyTreatmentData.of(null,
                        null,
                        readDrugsPostTreatment(treatmentData, treatmentCurator),
                        FormStatus.undefined()));
            }
        }
        return biopsyTreatmentDataList;
    }

    @NotNull
    private static List<DrugData> readDrugsPostTreatment(@NotNull WideTreatmentData treatmentData,
            @NotNull final TreatmentCurator treatmentCurator) {
        final List<DrugData> drugs = Lists.newArrayList();
        String drugName = treatmentData.drug();
        LocalDate drugsStartDate = treatmentData.startDate().isEmpty() ? null : createInterpretDate(treatmentData.startDate());
        LocalDate drugsEndDate = treatmentData.endDate().isEmpty() ? null : createInterpretDate(treatmentData.endDate());

        if (drugName != null || drugsStartDate != null || drugsEndDate != null) {
            final List<CuratedDrug> curatedDrugs = drugName == null ? Lists.newArrayList() : treatmentCurator.search(drugName);
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
                        createInterpretDate(responseData.date().isEmpty() ? null : responseData.date()),
                        null,
                        responseData.responseAccordingRecist(),
                        null,
                        FormStatus.undefined()));
            }
        }
        return biopsyTreatmentResponseDataList;
    }
}