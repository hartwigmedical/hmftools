package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeFormatterBuilder;
import java.util.List;
import java.util.Locale;

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
import com.hartwig.hmftools.patientdb.readers.wide.WideBiopsyData;
import com.hartwig.hmftools.patientdb.readers.wide.WideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.WidePreTreatmentData;
import com.hartwig.hmftools.patientdb.readers.wide.WideResponseData;
import com.hartwig.hmftools.patientdb.readers.wide.WideTreatmentData;

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
        // TODO: Create the timeline based on the wideEcrfModel
        return new Patient(patientIdentifier,
                toBaselineData(tumorLocationCurator.search(primaryTumorLocation)),
                preTreatmentData(wideEcrfModel.preTreatments()),
                sequencedSamples,
                toBiopsyData(wideEcrfModel.biopsies(), biopsySiteCurator),
                toBiopsyTreatmentData(wideEcrfModel.treatments(), treatmentCurator),
                toBiopsyTreatmentResponseData(wideEcrfModel.responses()),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList());
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
    private static PreTreatmentData preTreatmentData(@NotNull List<WidePreTreatmentData> widePreTreatmentData) {

        for (WidePreTreatmentData preTreatmentData : widePreTreatmentData) {

        }
        return ImmutablePreTreatmentData.of(null, null, Lists.newArrayList(), FormStatus.undefined());
    }

    @NotNull
    private static LocalDate createInterpretDate(@NotNull String date) {
        DateTimeFormatter inputFormatter =
                new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(Locale.ENGLISH);
        LocalDate localDate = LocalDate.parse(date, inputFormatter);

        DateTimeFormatter outputFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
        String formattedString = localDate.format(outputFormatter);
        return LocalDate.parse(formattedString);
    }

    @NotNull
    private static List<BiopsyData> toBiopsyData(@NotNull List<WideBiopsyData> wideBiopsyData,
            @NotNull BiopsySiteCurator biopsySiteCurator) {
        List<BiopsyData> biopsyDataList = Lists.newArrayList();
        final CuratedBiopsyType curatedBiopsyType = biopsySiteCurator.search(Strings.EMPTY, Strings.EMPTY, Strings.EMPTY, Strings.EMPTY);

        for (WideBiopsyData biopsyData : wideBiopsyData) {
            biopsyDataList.add(ImmutableBiopsyData.of(createInterpretDate(biopsyData.bioptDate()),
                    null,
                    null,
                    curatedBiopsyType,
                    null,
                    null,
                    FormStatus.undefined()));
        }
        return biopsyDataList;

    }

    @NotNull
    private static List<BiopsyTreatmentData> toBiopsyTreatmentData(@NotNull List<WideTreatmentData> wideTreatmentData,
            @NotNull final TreatmentCurator treatmentCurator) {
        List<BiopsyTreatmentData> biopsyTreatmentDataList = Lists.newArrayList();
        for (WideTreatmentData treatmentData : wideTreatmentData) {
            biopsyTreatmentDataList.add(ImmutableBiopsyTreatmentData.of(null,
                    null,
                    readDrugs(treatmentData, treatmentCurator),
                    FormStatus.undefined()));

        }
        return biopsyTreatmentDataList;
    }

    @NotNull
    public static List<DrugData> readDrugs(@NotNull WideTreatmentData treatmentData, @NotNull final TreatmentCurator treatmentCurator) {
        final List<DrugData> drugs = Lists.newArrayList();
        String drugName = treatmentData.drug();
        LocalDate drugsStartDate = createInterpretDate(treatmentData.startDate());
        LocalDate drugsEndDate = createInterpretDate(treatmentData.endDate());

        if (drugName != null || drugsStartDate != null || drugsEndDate != null) {
            final List<CuratedDrug> curatedDrugs = drugName == null ? Lists.newArrayList() : treatmentCurator.search(drugName);
            drugs.add(ImmutableDrugData.of(drugName, drugsStartDate, drugsEndDate, null, curatedDrugs));
        }
        return drugs;
    }

    @NotNull
    private static List<BiopsyTreatmentResponseData> toBiopsyTreatmentResponseData(@NotNull List<WideResponseData> wideResponseData) {
        List<BiopsyTreatmentResponseData> biopsyTreatmentResponseDataList = Lists.newArrayList();
        for (WideResponseData responseData : wideResponseData) {
            biopsyTreatmentResponseDataList.add(ImmutableBiopsyTreatmentResponseData.of(null,
                    createInterpretDate(responseData.date()),
                    null,
                    responseData.responseAccordingRecist(),
                    null,
                    FormStatus.undefined()));
        }
        return biopsyTreatmentResponseDataList;
    }
}