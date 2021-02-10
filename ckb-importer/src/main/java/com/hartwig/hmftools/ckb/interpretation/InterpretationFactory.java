package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodelinterpretation.CkbEntryInterpretation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.interpretation.clinicaltrial.ClinicalTrialFactory;
import com.hartwig.hmftools.ckb.interpretation.evidence.EvidenceFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class InterpretationFactory {

    private static final Logger LOGGER = LogManager.getLogger(InterpretationFactory.class);

    private InterpretationFactory() {
    }

    public static List<CkbEntryInterpretation> interpretationCkbDataModel(@NotNull CkbJsonDatabase ckbEntry) {
        List<CkbEntryInterpretation> CkbEntryInterpretation = Lists.newArrayList();
        int ckbId = 0;
        for (MolecularProfile molecularProfile : ckbEntry.molecularProfiles()) {
            ++ckbId;
            ImmutableCkbEntryInterpretation.Builder outputBuilder = ImmutableCkbEntryInterpretation.builder();
            outputBuilder.id(ckbId);

            ClinicalTrialFactory.interpretClinicalTrials(molecularProfile, ckbEntry, outputBuilder);

            EvidenceFactory.interpretEvidence(molecularProfile, ckbEntry, outputBuilder);

            LOGGER.info(outputBuilder.build());
            CkbEntryInterpretation.add(outputBuilder.build());
        }
        return CkbEntryInterpretation;
    }

}