package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.interpretation.clinicaltrial.ClinicalTrialFactory;
import com.hartwig.hmftools.ckb.interpretation.evidence.EvidenceFactory;
import com.hartwig.hmftools.ckb.interpretation.knowngenomicalteration.KnownGenomicAlterationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CkbEntryInterpretationFactory {

    private static final Logger LOGGER = LogManager.getLogger(CkbEntryInterpretationFactory.class);

    private CkbEntryInterpretationFactory() {
    }

    public static List<CkbEntryInterpretation> interpretationCkbDataModel(@NotNull CkbJsonDatabase ckbEntry) {
        List<CkbEntryInterpretation> CkbEntryInterpretation = Lists.newArrayList();
        int ckbId = 1;
        for (MolecularProfile molecularProfile : ckbEntry.molecularProfiles()) {
            ImmutableCkbEntryInterpretation.Builder outputBuilder = ImmutableCkbEntryInterpretation.builder();
            outputBuilder.id(ckbId);
            outputBuilder.molecularProfileId(molecularProfile.id());
            outputBuilder.profileName(molecularProfile.profileName());
            outputBuilder.createDate(molecularProfile.createDate());
            outputBuilder.updateDate(molecularProfile.updateDate());

            ClinicalTrialFactory.interpretClinicalTrials(molecularProfile, ckbEntry, outputBuilder);
            EvidenceFactory.interpretVariantEvidence(molecularProfile, ckbEntry, outputBuilder);
            KnownGenomicAlterationFactory.extractKnownGenomicAlteration(molecularProfile, ckbEntry, outputBuilder);

            LOGGER.info(outputBuilder.build());  //TODO removed when model is finished
            CkbEntryInterpretation.add(outputBuilder.build());
            ++ckbId;
        }
        return CkbEntryInterpretation;
    }
}