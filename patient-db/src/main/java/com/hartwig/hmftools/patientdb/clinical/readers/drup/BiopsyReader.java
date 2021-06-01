package com.hartwig.hmftools.patientdb.clinical.readers.drup;

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

import org.jetbrains.annotations.NotNull;

class BiopsyReader {

    private static final String STUDY_BIOPSY = "SE.BIOPBL";
    private static final String FORM_BIOPSY = "FRM.BIOPT";
    private static final String ITEMGROUP_BIOPSY = "GRP.BIOPT";
    private static final String FIELD_BIOPSY_TAKEN = "FLD.TBTAKEN";

    private static final String ITEMGROUP_TUMOR_BIOPSY = "GRP.TB";
    private static final String FIELD_BIOPSY_DATE = "FLD.TBDAT";
    private static final String FIELD_SITE = "FLD.TBSITE";
    private static final String FIELD_SITE_OTHER = "FLD.TBSITEOTH";
    private static final String FIELD_LOCATION = "FLD.TBLOC";

    @NotNull
    private final BiopsySiteCurator biopsySiteCurator;

    BiopsyReader(@NotNull final BiopsySiteCurator biopsySiteCurator) {
        this.biopsySiteCurator = biopsySiteCurator;
    }

    @NotNull
    List<BiopsyData> read(@NotNull EcrfPatient patient, @NotNull CuratedPrimaryTumor curatedPrimaryTumor) {
        List<BiopsyData> biopsies = Lists.newArrayList();
        for (EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BIOPSY)) {
            for (EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_BIOPSY)) {
                String biopsyTaken = null;
                // We assume 1:1 relation between biopsy form and tumor biopsy form
                for (EcrfItemGroup biopsyGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_BIOPSY)) {
                    biopsyTaken = biopsyGroup.readItemString(FIELD_BIOPSY_TAKEN);
                }

                for (EcrfItemGroup tumorBiopsyGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_TUMOR_BIOPSY)) {
                    LocalDate date = tumorBiopsyGroup.readItemDate(FIELD_BIOPSY_DATE);

                    String site = tumorBiopsyGroup.readItemString(FIELD_SITE);
                    String siteOther = tumorBiopsyGroup.readItemString(FIELD_SITE_OTHER);
                    String finalSite = (site == null || site.trim().toLowerCase().startsWith("other")) ? siteOther : site;

                    String location = tumorBiopsyGroup.readItemString(FIELD_LOCATION);

                    CuratedBiopsyType curatedBiopsyType = biopsySiteCurator.search(curatedPrimaryTumor.location(),
                            curatedPrimaryTumor.type(),
                            finalSite,
                            location);

                    biopsies.add(BiopsyData.of(date, biopsyTaken, null, curatedBiopsyType, finalSite, location, form.status()));
                }
            }
        }

        return biopsies;
    }
}
