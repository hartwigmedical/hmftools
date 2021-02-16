package com.hartwig.hmftools.ckb.datamodel;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialFactory;
import com.hartwig.hmftools.ckb.datamodel.evidence.EvidenceFactory;
import com.hartwig.hmftools.ckb.datamodel.variant.VariantFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.molecularprofile.JsonMolecularProfile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class JsonDatabaseToCkbEntryConverter {

    private static final Logger LOGGER = LogManager.getLogger(JsonDatabaseToCkbEntryConverter.class);

    private JsonDatabaseToCkbEntryConverter() {
    }

    @NotNull
    public static List<CkbEntry> convert(@NotNull CkbJsonDatabase ckbJsonDatabase) {
        List<CkbEntry> ckbEntries = Lists.newArrayList();
        for (JsonMolecularProfile molecularProfile : ckbJsonDatabase.molecularProfiles()) {
            ckbEntries.add(ImmutableCkbEntry.builder()
                    .profileId(molecularProfile.id())
                    .profileName(molecularProfile.profileName())
                    .createDate(molecularProfile.createDate())
                    .updateDate(molecularProfile.updateDate())
                    .variants(VariantFactory.extractVariants(ckbJsonDatabase, molecularProfile.geneVariants()))
                    .evidences(EvidenceFactory.extractEvidences(ckbJsonDatabase,
                            molecularProfile.variantLevelEvidence().evidences(),
                            molecularProfile.id()))
                    .clinicalTrials(ClinicalTrialFactory.extractClinicalTrials(ckbJsonDatabase,
                            molecularProfile.variantAssociatedClinicalTrials()))
                    .build());
        }

        return ckbEntries;
    }
}