package com.hartwig.hmftools.ckb.interpretation.evidence;

import com.hartwig.hmftools.ckb.datamodelinterpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;

import org.jetbrains.annotations.NotNull;

public class EvidenceFactory {

    private EvidenceFactory() {

    }

    public static void interpretEvidence(@NotNull MolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntryInterpretation.Builder outputBuilder) {
        for (EvidenceInfo evidenceInfo : molecularProfile.variantLevelEvidence().evidence()) {
            outputBuilder.addEvidenceInterpretation();
        }
    }
}
