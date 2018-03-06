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
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.CuratedTreatment;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.ImmutableDrugData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BiopsyTreatmentReader {
    private static final String STUDY_AFTERBIOPT = "SE.AFTERBIOPT";
    public static final String FORM_TREATMENT = "FRM.TRTAFTER";
    private static final String ITEMGROUP_TREATMENT_AFTER = "GRP.TRTAFTER.TRTAFTER";
    private static final String ITEMGROUP_SYSPOSTBIO = "GRP.TRTAFTER.SYSPOSTBIO";
    public static final String FIELD_TREATMENT_GIVEN = "FLD.TRTAFTER.SYSTEMICST";
    public static final String FIELD_DRUG_START = "FLD.TRTAFTER.SYSSTDT";
    public static final String FIELD_DRUG_END = "FLD.TRTAFTER.SYSENDT";
    public static final String FIELD_DRUG = "FLD.TRTAFTER.PLANNEDTRT";
    public static final String FIELD_DRUG_OTHER = "FLD.TRTAFTER.SYSREGPOST";

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final TreatmentCurator treatmentCurator;

    BiopsyTreatmentReader(@NotNull final TreatmentCurator treatmentCurator) {
        this.treatmentCurator = treatmentCurator;
    }

    @NotNull
    List<BiopsyTreatmentData> read(@NotNull final EcrfPatient patient) throws IOException {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_AFTERBIOPT)) {
            for (final EcrfForm treatmentForm : studyEvent.nonEmptyFormsPerOID(FORM_TREATMENT, false)) {
                final String treatmentGiven = readTreatmentGiven(treatmentForm);
                final List<DrugData> drugs = readDrugs(treatmentForm);
                treatments.add(ImmutableBiopsyTreatmentData.of(treatmentGiven, drugs, treatmentForm.status(), treatmentForm.locked()));
            }
        }
        return treatments;
    }

    @NotNull
    private List<DrugData> readDrugs(@NotNull final EcrfForm treatmentForm) throws IOException {
        final List<DrugData> drugs = Lists.newArrayList();
        for (final EcrfItemGroup itemGroup : treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_SYSPOSTBIO, false)) {
            final LocalDate drugStart = itemGroup.readItemDate(FIELD_DRUG_START, 0, DATE_FORMATTER, false);
            final LocalDate drugEnd = itemGroup.readItemDate(FIELD_DRUG_END, 0, DATE_FORMATTER, false);
            String drugName = itemGroup.readItemString(FIELD_DRUG, 0, false);
            if (drugName == null || drugName.trim().toLowerCase().startsWith("other")) {
                drugName = itemGroup.readItemString(FIELD_DRUG_OTHER, 0, false);
            }
            final List<CuratedTreatment> curatedDrugs = drugName == null ? Lists.newArrayList() : treatmentCurator.search(drugName);
            drugs.add(ImmutableDrugData.of(drugName, drugStart, drugEnd, null, curatedDrugs));
        }
        return drugs;
    }

    @Nullable
    private static String readTreatmentGiven(@NotNull final EcrfForm treatmentForm) {
        final List<EcrfItemGroup> itemGroups = treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_TREATMENT_AFTER, false);
        if (itemGroups.size() > 0) {
            return itemGroups.get(0).readItemString(FIELD_TREATMENT_GIVEN, 0, false);
        }
        return null;
    }
}
