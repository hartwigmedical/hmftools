package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;
import com.hartwig.hmftools.patientdb.data.CuratedTreatment;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentDrugData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BiopsyTreatmentReader {
    private static final Logger LOGGER = LogManager.getLogger(BiopsyTreatmentReader.class);

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
        final List<BiopsyTreatmentData> treatmentDatas = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_AFTERBIOPT)) {
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_TREATMENT, false)) {
                final String treatmentGiven = readTreatmentGiven(form);
                final List<BiopsyTreatmentDrugData> drugs = readDrugs(form);
                final LocalDate treatmentStart = determineTreatmentStartDate(drugs);
                final LocalDate treatmentEnd = determineTreatmentEndDate(drugs);
                treatmentDatas.add(
                        ImmutableBiopsyTreatmentData.of(treatmentGiven, treatmentStart, treatmentEnd, drugs, form.status(), form.locked()));
            }
        }
        return treatmentDatas;
    }

    @NotNull
    private List<BiopsyTreatmentDrugData> readDrugs(@NotNull final EcrfForm form) throws IOException {
        final List<BiopsyTreatmentDrugData> drugs = Lists.newArrayList();
        for (final EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_SYSPOSTBIO, false)) {
            final LocalDate drugStart = itemGroup.readItemDate(FIELD_DRUG_START, 0, DATE_FORMATTER, false);
            final LocalDate drugEnd = itemGroup.readItemDate(FIELD_DRUG_END, 0, DATE_FORMATTER, false);
            String drugName = itemGroup.readItemString(FIELD_DRUG, 0, false);
            if (drugName == null || drugName.trim().toLowerCase().startsWith("other")) {
                drugName = itemGroup.readItemString(FIELD_DRUG_OTHER, 0, false);
            }
            if (drugName == null) {
                drugs.add(ImmutableBiopsyTreatmentDrugData.of(null, null, drugStart, drugEnd));
            } else {
                drugs.addAll(curatedDrugs(drugName, drugStart, drugEnd));
            }
        }
        return drugs;
    }

    @NotNull
    private List<BiopsyTreatmentDrugData> curatedDrugs(@NotNull final String drugName, @Nullable final LocalDate drugStart,
            @Nullable final LocalDate drugEnd) throws IOException {
        final List<CuratedTreatment> matchedTreatments = treatmentCurator.search(drugName);
        if (matchedTreatments.isEmpty()) {
            LOGGER.warn(
                    "Failed to map ecrf drug {} to any drug in curated list. Either the list contained no entry, or search results were ambiguous.",
                    drugName);
            return Lists.newArrayList(ImmutableBiopsyTreatmentDrugData.of(drugName, null, drugStart, drugEnd));
        } else {
            return matchedTreatments.stream()
                    .map(curatedTreatment -> ImmutableBiopsyTreatmentDrugData.of(curatedTreatment.name(), curatedTreatment.type(),
                            drugStart, drugEnd))
                    .collect(Collectors.toList());
        }
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
        final List<EcrfItemGroup> itemGroups = form.nonEmptyItemGroupsPerOID(ITEMGROUP_TREATMENT_AFTER, false);
        if (itemGroups.size() > 0) {
            return itemGroups.get(0).readItemString(FIELD_TREATMENT_GIVEN, 0, false);
        }
        return null;
    }
}
