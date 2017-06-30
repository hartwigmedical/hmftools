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
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_BIOPS, true)) {
                for (final EcrfItemGroup itemGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_BIOPSIES, true)) {
                    final LocalDate date = itemGroup.readItemDate(FIELD_BIOPSY_DATE, 0, DATE_FORMATTER, true);
                    final String location = itemGroup.readItemString(FIELD_LOCATION, 0, true);
                    if (location != null && location.trim().toLowerCase().startsWith("other")) {
                        final String location_other = itemGroup.readItemString(FIELD_LOCATION_OTHER, 0, true);
                        biopsies.add(new BiopsyData(date, location_other));
                    } else {
                        biopsies.add(new BiopsyData(date, location));
                    }
                }
            }
        }
        return biopsies;
    }
}
