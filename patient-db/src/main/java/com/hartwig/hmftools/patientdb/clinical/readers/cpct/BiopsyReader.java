package com.hartwig.hmftools.patientdb.clinical.readers.cpct;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedBiopsyType;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfStudyEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class BiopsyReader {

    static final String STUDY_BIOPSY = "SE.BIOPSY";
    static final String FORM_BIOPS = "FRM.BIOPS";
    static final String ITEMGROUP_BIOPSY = "GRP.BIOPS.BIOPS";
    static final String FIELD_BIOPSY_TAKEN = "FLD.BIOPS.CPCT";
    static final String FIELD_BIOPSY_EVALUABLE = "FLD.BIOPS.BIOPEFS";

    static final String ITEMGROUP_BIOPSIES = "GRP.BIOPS.BIOPSIES";
    static final String FIELD_BIOPSY_DATE = "FLD.BIOPS.BIOPTDT";
    static final String FIELD_SITE = "FLD.BIOPS.BILESSITE";
    static final String FIELD_SITE_OTHER = "FLD.BIOPS.BIOTHLESSITE";
    static final String FIELD_LOCATION = "FLD.BIOPS.BILESLOC";

    @NotNull
    private final BiopsySiteCurator biopsySiteCurator;

    BiopsyReader(@NotNull final BiopsySiteCurator biopsySiteCurator) {
        this.biopsySiteCurator = biopsySiteCurator;
    }

    @NotNull
    List<BiopsyData> read(@NotNull EcrfPatient patient, @NotNull CuratedPrimaryTumor curatedPrimaryTumor) {
        List<BiopsyData> biopsies = Lists.newArrayList();
        for (EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BIOPSY)) {
            for (EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_BIOPS)) {
                String biopsyTaken = null;
                String biopsyEvaluable = null;
                // This works as there is generally a 1:N relation between BIOPSY and BIOPSIES item groups.
                for (EcrfItemGroup biopsyGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_BIOPSY)) {
                    biopsyTaken = biopsyGroup.readItemString(FIELD_BIOPSY_TAKEN);
                    biopsyEvaluable = biopsyGroup.readItemString(FIELD_BIOPSY_EVALUABLE);
                }

                for (EcrfItemGroup biopsiesGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_BIOPSIES)) {
                    LocalDate date = biopsiesGroup.readItemDate(FIELD_BIOPSY_DATE);

                    String site = biopsiesGroup.readItemString(FIELD_SITE);
                    String siteOther = biopsiesGroup.readItemString(FIELD_SITE_OTHER);
                    String finalSite = (site == null || site.trim().toLowerCase().startsWith("other")) ? siteOther : site;

                    String location = biopsiesGroup.readItemString(FIELD_LOCATION);

                    String subTumorType = Strings.EMPTY;
                    if (curatedPrimaryTumor.type() != null) {
                        subTumorType = curatedPrimaryTumor.subType().equals(Strings.EMPTY)
                                ? curatedPrimaryTumor.type()
                                : curatedPrimaryTumor.subType();
                    }

                    CuratedBiopsyType curatedBiopsyType = biopsySiteCurator.search(
                            curatedPrimaryTumor.location(),
                            subTumorType,
                            finalSite,
                            location);

                    BiopsyData biopsy =
                            BiopsyData.of(date, biopsyTaken, biopsyEvaluable, curatedBiopsyType, finalSite, location, form.status());

                    // The ecrf contains many duplicate forms that are impossible to remove. This is because in the past a new biopsy
                    // form needed to be created for every treatment response.
                    if (!isDuplicate(biopsies, biopsy) && !isEmpty(biopsy)) {
                        biopsies.add(biopsy);
                    }
                }
            }
        } return biopsies;
    }

    private static boolean isEmpty(@NotNull BiopsyData biopsy) {
        return (biopsy.date() == null && biopsy.location() == null && biopsy.site() == null);
    }

    private static boolean isDuplicate(@NotNull List<BiopsyData> biopsies, @NotNull BiopsyData biopsyToCheck) {
        for (BiopsyData biopsy : biopsies) {
            if (equals(biopsy.date(), biopsyToCheck.date()) && equals(biopsy.location(), biopsyToCheck.location()) && equals(biopsy.site(),
                    biopsyToCheck.site())) {
                return true;
            }
        }
        return false;
    }

    private static boolean equals(@Nullable Object field1, @Nullable Object field2) {
        return ((field1 == null && field2 == null) || (field1 != null && field1.equals(field2)));
    }
}
