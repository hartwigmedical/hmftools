package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyData;

import org.jetbrains.annotations.NotNull;

public final class BiopsyReader {
    private static final String STUDY_BIOPSY = "SE.BIOPSY";
    public static final String FORM_BIOPS = "FRM.BIOPS";
    private static final String ITEMGROUP_BIOPSIES = "GRP.BIOPS.BIOPSIES";

    public static final String FIELD_BIOPSY_DATE = "FLD.BIOPS.BIOPTDT";
    public static final String FIELD_LOCATION = "FLD.BIOPS.BILESSITE";
    public static final String FIELD_LOCATION_OTHER = "FLD.BIOPS.BIOTHLESSITE";

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    private BiopsyReader() {
    }

    @NotNull
    static List<BiopsyData> read(@NotNull final EcrfPatient patient) {
        final List<BiopsyData> biopsies = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BIOPSY)) {
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_BIOPS, false)) {
                for (final EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_BIOPSIES, false)) {
                    final LocalDate date = itemGroup.readItemDate(FIELD_BIOPSY_DATE, 0, DATE_FORMATTER, false);
                    final String location = itemGroup.readItemString(FIELD_LOCATION, 0, false);
                    if (location == null || location.trim().toLowerCase().startsWith("other")) {
                        final String location_other = itemGroup.readItemString(FIELD_LOCATION_OTHER, 0, false);
                        biopsies.add(ImmutableBiopsyData.of(date, location_other, form.status(), form.locked()));
                    } else {
                        biopsies.add(ImmutableBiopsyData.of(date, location, form.status(), form.locked()));
                    }
                }
            }
        }
        return biopsies;
    }
}
