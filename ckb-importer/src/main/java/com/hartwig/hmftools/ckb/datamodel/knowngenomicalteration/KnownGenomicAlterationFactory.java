package com.hartwig.hmftools.ckb.datamodel.knowngenomicalteration;

import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.common.molecularprofile.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.molecularprofile.JsonMolecularProfile;

import org.jetbrains.annotations.NotNull;

public final class KnownGenomicAlterationFactory {

    private KnownGenomicAlterationFactory() {
    }

    public static void extractKnownGenomicAlteration(@NotNull JsonMolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntry.Builder outputBuilder) {
        outputBuilder.variants(MolecularProfileInterpretationFactory.extractVariantGeneInfo(ckbEntry, molecularProfile));
    }
}