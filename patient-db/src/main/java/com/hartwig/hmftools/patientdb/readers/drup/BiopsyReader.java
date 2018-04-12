package com.hartwig.hmftools.patientdb.readers.drup;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.CuratedBiopsyType;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyData;

import org.jetbrains.annotations.NotNull;

class BiopsyReader {

    private static final String STUDY_BIOPSY = "SE.BIOPBL";
    private static final String FORM_BIOPSY = "FRM.BIOPT";
    private static final String ITEMGROUP_BIOPSY = "GRP.BIOPT.BIOPT";
    private static final String FIELD_BIOPSY_TAKEN = "FLD.BIOPT_TBTAKEN";

    private static final String ITEMGROUP_TUMOR_BIOPSY = "GRP.BIOPT.TB";
    private static final String FIELD_BIOPSY_DATE = "FLD.BIOPT.TBDAT";
    private static final String FIELD_SITE = "FLD.BIOPT.TBSITE";
    private static final String FIELD_SITE_OTHER = "FLD.BIOPT.TBSITEOTH";
    private static final String FIELD_LOCATION = "FLD.BIOPT.TBLOC";

    @NotNull
    private final BiopsySiteCurator biopsySiteCurator;

    BiopsyReader(@NotNull final BiopsySiteCurator biopsySiteCurator) {
        this.biopsySiteCurator = biopsySiteCurator;
    }

    @NotNull
    List<BiopsyData> read(@NotNull EcrfPatient patient, @NotNull CuratedTumorLocation curatedTumorLocation) {
        final List<BiopsyData> biopsies = Lists.newArrayList();
        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BIOPSY)) {
            for (final EcrfForm form : studyEvent.nonEmptyFormsPerOID(FORM_BIOPSY)) {
                String biopsyTaken = null;
                // KODU: We assume 1:1 relation between biopsy form and tumor biopsy form
                for (final EcrfItemGroup biopsyGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_BIOPSY)) {
                    biopsyTaken = biopsyGroup.readItemString(FIELD_BIOPSY_TAKEN);
                }

                for (final EcrfItemGroup tumorBiopsyGroup : form.nonEmptyItemGroupsPerOID(ITEMGROUP_TUMOR_BIOPSY)) {
                    final LocalDate date = tumorBiopsyGroup.readItemDate(FIELD_BIOPSY_DATE);

                    final String site = tumorBiopsyGroup.readItemString(FIELD_SITE);
                    final String siteOther = tumorBiopsyGroup.readItemString(FIELD_SITE_OTHER);
                    final String finalSite = (site == null || site.trim().toLowerCase().startsWith("other")) ? siteOther : site;

                    final String location = tumorBiopsyGroup.readItemString(FIELD_LOCATION);

                    final CuratedBiopsyType curatedBiopsyType = biopsySiteCurator.search(curatedTumorLocation.primaryTumorLocation(),
                            curatedTumorLocation.subType(),
                            finalSite,
                            location);

                    biopsies.add(ImmutableBiopsyData.of(date, biopsyTaken, null, curatedBiopsyType, finalSite, location, form.status()));
                }
            }
        }

        return biopsies;
    }
}
