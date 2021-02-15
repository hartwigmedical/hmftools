package com.hartwig.hmftools.ckb.datamodel;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialFactory;
import com.hartwig.hmftools.ckb.datamodel.evidence.EvidenceFactory;
import com.hartwig.hmftools.ckb.datamodel.knowngenomicalteration.KnownGenomicAlterationFactory;
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
            ImmutableCkbEntry.Builder outputBuilder = ImmutableCkbEntry.builder();
            outputBuilder.profileId(molecularProfile.id());
            outputBuilder.profileName(molecularProfile.profileName());
            outputBuilder.createDate(molecularProfile.createDate());
            outputBuilder.updateDate(molecularProfile.updateDate());

            ClinicalTrialFactory.interpretClinicalTrials(molecularProfile, ckbJsonDatabase, outputBuilder);
            EvidenceFactory.interpretVariantEvidence(molecularProfile, ckbJsonDatabase, outputBuilder);
            KnownGenomicAlterationFactory.extractKnownGenomicAlteration(molecularProfile, ckbJsonDatabase, outputBuilder);

            LOGGER.info(outputBuilder.build());  //TODO removed when model is finished
            ckbEntries.add(outputBuilder.build());
        }

        return ckbEntries;
    }
}