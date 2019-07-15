package com.hartwig.hmftools.patientreporter.genepanel;

import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;

import org.jetbrains.annotations.NotNull;

public final class GeneModelFactory {

    private GeneModelFactory() {
    }

    @NotNull
    public static GeneModel create() {
        return ImmutableGeneModel.of(DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet(),
                DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet(),
                CNADrivers.amplificationTargets(),
                CNADrivers.deletionTargets());
    }
}
