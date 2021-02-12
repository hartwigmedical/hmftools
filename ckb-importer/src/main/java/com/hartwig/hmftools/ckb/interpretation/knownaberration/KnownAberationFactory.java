package com.hartwig.hmftools.ckb.interpretation.knownaberration;

import com.hartwig.hmftools.ckb.datamodelinterpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.variantinterpretation.VariantInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;

import org.jetbrains.annotations.NotNull;

public class KnownAberationFactory {

    private KnownAberationFactory() {
    }

    public static void extractKnownGenomicAberations(@NotNull MolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntryInterpretation.Builder outputBuilder) {

        outputBuilder.knownAberations(ImmutableKnownAberation.builder()
                .knownAberations(VariantInterpretationFactory.extractVariantGeneInfo(ckbEntry, molecularProfile).build())
                .build());
    }
}
