package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.data.CuratedTreatment;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableDrugData;
import com.hartwig.hmftools.patientdb.data.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PreTreatmentReader {

    private static final Logger LOGGER = LogManager.getLogger(PreTreatmentReader.class);

    private static final String STUDY_BASELINE = "SE.BASELINE";

    private static final String FORM_TREATMENT = "FRM.PRETHERAPY";
    private static final String ITEMGROUP_TREATMENT = "GRP.PRETHERAPY.PRETHERAPY";
    private static final String FIELD_TREATMENT_GIVEN = "FLD.PRETHERAPY.SYSTEMIC";
    private static final String FIELD_RADIOTHERAPY_GIVEN = "FLD.PRETHERAPY.RADIOTHER";
    private static final String ITEMGROUP_DRUGS = "GRP.PRETHERAPY.SYSTEMICTRT";
    private static final String FIELD_PRE_DRUG_START = "FLD.PRETHERAPY.SYSTEMICSTDTC";
    private static final String FIELD_PRE_DRUG_END = "FLD.PRETHERAPY.SYSTEMICENDTC";
    public static final String FIELD_PRE_DRUG = "FLD.PRETHERAPY.SYSTEMICREG";
    private static final String FIELD_PRE_BEST_RESPONSE = "FLD.PRETHERAPY.SYSTEMICRESP";

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final TreatmentCurator treatmentCurator;

    PreTreatmentReader(@NotNull final TreatmentCurator treatmentCurator) {
        this.treatmentCurator = treatmentCurator;
    }

    @NotNull
    PreTreatmentData read(@NotNull final EcrfPatient patient) throws IOException {
        PreTreatmentData preTreatmentData = null;
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BASELINE)) {
            for (final EcrfForm treatmentForm : studyEvent.nonEmptyFormsPerOID(FORM_TREATMENT, false)) {
                String treatmentGiven = readTreatmentGiven(treatmentForm);
                String radiotherapyGiven = readRadiotherapyGiven(treatmentForm);
                final List<DrugData> drugs = readDrugs(treatmentForm);
                if (preTreatmentData == null) {
                    preTreatmentData = ImmutablePreTreatmentData.of(treatmentGiven,
                            radiotherapyGiven,
                            drugs,
                            treatmentForm.status(),
                            treatmentForm.locked());
                } else {
                    LOGGER.warn("Multiple pre-therapy forms for found patient: " + patient.patientId());
                }
            }
        }

        return preTreatmentData != null
                ? preTreatmentData
                : ImmutablePreTreatmentData.of(null, null, Lists.newArrayList(), FormStatusState.UNKNOWN, false);
    }

    @NotNull
    private List<DrugData> readDrugs(@NotNull final EcrfForm treatmentForm) throws IOException {
        final List<DrugData> drugs = Lists.newArrayList();
        for (final EcrfItemGroup itemGroup : treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_DRUGS, false)) {
            final LocalDate drugStart = itemGroup.readItemDate(FIELD_PRE_DRUG_START, 0, DATE_FORMATTER, false);
            final LocalDate drugEnd = itemGroup.readItemDate(FIELD_PRE_DRUG_END, 0, DATE_FORMATTER, false);
            final String drugName = itemGroup.readItemString(FIELD_PRE_DRUG, 0, false);
            final String bestResponse = itemGroup.readItemString(FIELD_PRE_BEST_RESPONSE, 0, false);

            final List<CuratedTreatment> curatedDrugs = drugName == null ? Lists.newArrayList() : treatmentCurator.search(drugName);
            drugs.add(ImmutableDrugData.of(drugName, drugStart, drugEnd, bestResponse, curatedDrugs));
        }
        return drugs;
    }

    @Nullable
    private static String readTreatmentGiven(@NotNull final EcrfForm treatmentForm) {
        final List<EcrfItemGroup> itemGroups = treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_TREATMENT, false);
        if (itemGroups.size() > 0) {
            return itemGroups.get(0).readItemString(FIELD_TREATMENT_GIVEN, 0, false);
        }
        return null;
    }

    @Nullable
    private static String readRadiotherapyGiven(@NotNull final EcrfForm treatmentForm) {
        final List<EcrfItemGroup> itemGroups = treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_TREATMENT, false);
        if (itemGroups.size() > 0) {
            return itemGroups.get(0).readItemString(FIELD_RADIOTHERAPY_GIVEN, 0, false);
        }
        return null;
    }
}
