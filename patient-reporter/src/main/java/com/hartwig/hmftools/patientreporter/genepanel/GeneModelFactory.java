package com.hartwig.hmftools.patientreporter.genepanel;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.drup.DrupActionabilityModel;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.jetbrains.annotations.NotNull;

public final class GeneModelFactory {

    private static final Set<String> DISRUPTION_GENE_WHITE_LIST = Sets.newHashSet("ALK", "NTRK1", "NTRK2", "NTRK3", "RET", "ROS1", "NRG1");

    private GeneModelFactory() {
    }

    @NotNull
    public static GeneModel create(@NotNull DrupActionabilityModel drupActionabilityModel) {
        Map<String, DriverCategory> geneDriverCategoryMap = Maps.newHashMap();

        Set<String> somaticVariantDriverGenes = Sets.newHashSet();

        for (String oncoGene : DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet()) {
            somaticVariantDriverGenes.add(oncoGene);
            geneDriverCategoryMap.put(oncoGene, DriverCategory.ONCO);
        }

        for (String tsgGene : DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet()) {
            somaticVariantDriverGenes.add(tsgGene);
            geneDriverCategoryMap.put(tsgGene, DriverCategory.TSG);
        }

        // For any gene that is actionable in DRUP and configured to onco/tsg, use this in driver category map,
        // only if not configured in the driver catalog.
        for (Map.Entry<String, DriverCategory> drupDriverCategoryEntry : drupActionabilityModel.geneDriverCategoryMap().entrySet()) {
            geneDriverCategoryMap.putIfAbsent(drupDriverCategoryEntry.getKey(), drupDriverCategoryEntry.getValue());
        }

        return ImmutableGeneModel.of(somaticVariantDriverGenes,
                CNADrivers.amplificationTargets(),
                CNADrivers.deletionTargets(),
                drupActionabilityModel.actionableGenes(),
                DISRUPTION_GENE_WHITE_LIST,
                geneDriverCategoryMap);
    }
}
