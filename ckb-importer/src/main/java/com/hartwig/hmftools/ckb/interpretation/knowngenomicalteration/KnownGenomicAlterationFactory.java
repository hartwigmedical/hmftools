package com.hartwig.hmftools.ckb.interpretation.knowngenomicalteration;

import com.hartwig.hmftools.ckb.interpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.molecularprofileinterpretation.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;

import org.jetbrains.annotations.NotNull;

public class KnownGenomicAlterationFactory {

    private KnownGenomicAlterationFactory() {
    }

    public static void extractKnownGenomicAberations(@NotNull MolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntryInterpretation.Builder outputBuilder) {

        outputBuilder.knownAberation(ImmutableKnownGenomicAlteration.builder()
                .knownAberations(MolecularProfileInterpretationFactory.extractVariantGeneInfo(ckbEntry, molecularProfile).build())
                .build());
    }
}