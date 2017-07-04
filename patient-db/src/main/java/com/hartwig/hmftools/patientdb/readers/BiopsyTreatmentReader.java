package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;

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
    private final Map<String, String> treatmentToTypeMapping;

    BiopsyTreatmentReader(@NotNull final Map<String, String> treatmentToTypeMapping) {
        this.treatmentToTypeMapping = treatmentToTypeMapping;
    }

    @NotNull
    List<BiopsyTreatmentData> read(@NotNull final EcrfPatient patient) {
        final List<BiopsyTreatmentData> treatmentDatas = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_AFTERBIOPT)) {
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_TREATMENT, true)) {
                final String treatmentGiven = readTreatmentGiven(form);
                final List<BiopsyTreatmentDrugData> drugs = readDrugs(form);
                final LocalDate treatmentStart = determineTreatmentStartDate(drugs);
                final LocalDate treatmentEnd = determineTreatmentEndDate(drugs);
                treatmentDatas.add(new BiopsyTreatmentData(treatmentGiven, treatmentStart, treatmentEnd, drugs));
            }
        }
        return treatmentDatas;
    }

    @NotNull
    private List<BiopsyTreatmentDrugData> readDrugs(@NotNull final EcrfForm form) {
        final List<BiopsyTreatmentDrugData> drugs = Lists.newArrayList();
        for (final EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_SYSPOSTBIO, true)) {
            final LocalDate drugStart = itemGroup.readItemDate(FIELD_DRUG_START, 0, DATE_FORMATTER, true);
            final LocalDate drugEnd = itemGroup.readItemDate(FIELD_DRUG_END, 0, DATE_FORMATTER, true);
            String drugName = itemGroup.readItemString(FIELD_DRUG, 0, true);
            if (drugName == null || drugName.trim().toLowerCase().startsWith("other")) {
                drugName = itemGroup.readItemString(FIELD_DRUG_OTHER, 0, true);
            }
            final String drugType = drugName == null ? null : treatmentToTypeMapping.get(drugName.toLowerCase().trim());
            drugs.add(new BiopsyTreatmentDrugData(drugName, drugType, drugStart, drugEnd));
        }
        return drugs;
    }

    @Nullable
    private static LocalDate determineTreatmentStartDate(@NotNull final List<BiopsyTreatmentDrugData> drugs) {
        LocalDate startDate = null;
        for (final BiopsyTreatmentDrugData drug : drugs) {
            final LocalDate drugStartDate = drug.startDate();
            if (startDate == null || (drugStartDate != null && drugStartDate.isBefore(startDate))) {
                startDate = drugStartDate;
            }
        }
        return startDate;
    }

    @Nullable
    private static LocalDate determineTreatmentEndDate(@NotNull final List<BiopsyTreatmentDrugData> drugs) {
        if (drugs.isEmpty()) {
            return null;
        } else {
            LocalDate endDate = drugs.get(0).endDate();
            for (final BiopsyTreatmentDrugData drug : drugs) {
                final LocalDate drugEndDate = drug.endDate();
                if (drugEndDate == null || (endDate != null && drugEndDate.isAfter(endDate))) {
                    endDate = drugEndDate;
                }
            }
            return endDate;
        }
    }

    @Nullable
    private static String readTreatmentGiven(@NotNull final EcrfForm form) {
        final List<EcrfItemGroup> itemGroups = form.nonEmptyItemGroupsPerOID(ITEMGROUP_TREATMENT_AFTER, true);
        if (itemGroups.size() > 0) {
            return itemGroups.get(0).readItemString(FIELD_TREATMENT_GIVEN, 0, true);
        }
        return null;
    }
}
