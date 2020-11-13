package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;

class DndsGeneName {

    private final DndsGeneNameMap dndsGeneNameMap;
    private final Predicate<DriverGene> alignmentMapPredicate;
    private final Function<DriverGene, String> dndsGeneFunction;
    private final Set<String> dndsGenes;

    public DndsGeneName(final DriverGenePanelAssembly assembly, final Set<String> dndsGenes) {
        this.dndsGenes = dndsGenes;
        dndsGeneNameMap = new DndsGeneNameMap();
        if (assembly == DriverGenePanelAssembly.HG19) {
            alignmentMapPredicate = aDriverGene -> dndsGeneNameMap.isValidHg19Gene(aDriverGene.gene());
            dndsGeneFunction = DriverGene::gene;
        } else {
            alignmentMapPredicate = aDriverGene -> dndsGeneNameMap.isValidHg38Gene(aDriverGene.gene());
            dndsGeneFunction = aDriverGene -> dndsGeneNameMap.hg19Gene(aDriverGene.gene());
        }
    }

    public boolean isValid(DriverGene driverGene) {
        if (!alignmentMapPredicate.test(driverGene)) {
            return false;
        }

        String dndsGeneName = dndsGeneName(driverGene);
        return dndsGenes.contains(dndsGeneName);
    }

    public String dndsGeneName(DriverGene driverGene) {
        return dndsGeneFunction.apply(driverGene);
    }
}
