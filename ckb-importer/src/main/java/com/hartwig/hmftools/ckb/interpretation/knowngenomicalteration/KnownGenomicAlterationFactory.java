package com.hartwig.hmftools.ckb.interpretation.knowngenomicalteration;

import com.hartwig.hmftools.ckb.interpretation.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.interpretation.common.molecularprofileinterpretation.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;

import org.jetbrains.annotations.NotNull;

public class KnownGenomicAlterationFactory {

    private KnownGenomicAlterationFactory() {
    }

    public static void extractKnownGenomicAlteration(@NotNull MolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntry.Builder outputBuilder) {

        outputBuilder.knownGenomicAlteration(ImmutableKnownGenomicAlteration.builder()
                .knownGenomicAlterationInterpretation(MolecularProfileInterpretationFactory.extractVariantGeneInfo(ckbEntry,
                        molecularProfile).build())
                .build());
    }
}