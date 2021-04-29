package com.hartwig.hmftools.serve;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;

import org.jetbrains.annotations.NotNull;

public final class DriverGeneTestFactory {

    private DriverGeneTestFactory() {
    }

    @NotNull
    public static List<DriverGene> createDriverGenes(@NotNull String geneTsg, @NotNull String geneOnco) {
        ImmutableDriverGene.Builder driverGeneBuilder = ImmutableDriverGene.builder()
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(false)
                .reportAmplification(false)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE);

        DriverGene driverGeneTsg = driverGeneBuilder.gene(geneTsg).likelihoodType(DriverCategory.TSG).build();
        DriverGene driverGeneOnco = driverGeneBuilder.gene(geneOnco).likelihoodType(DriverCategory.ONCO).build();

        return Lists.newArrayList(driverGeneTsg, driverGeneOnco);
    }
}
