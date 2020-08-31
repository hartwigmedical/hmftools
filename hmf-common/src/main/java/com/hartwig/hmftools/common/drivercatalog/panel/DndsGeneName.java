package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.function.Function;
import java.util.function.Predicate;

class DndsGeneName {

    private final DriverGenePanelMap panelMap;
    private final Predicate<DriverGene> isValidPredicate;
    private final Function<DriverGene, String> dndsGeneFunction;

    public DndsGeneName(final DriverGenePanelAssembly version) {
        panelMap = new DriverGenePanelMap();
        if (version == DriverGenePanelAssembly.HG19) {
            isValidPredicate = aDriverGene -> panelMap.isValidHg19Gene(aDriverGene.gene());
            dndsGeneFunction = DriverGene::gene;
        } else {
            isValidPredicate = aDriverGene -> panelMap.isValidHg38Gene(aDriverGene.gene());
            dndsGeneFunction = aDriverGene -> panelMap.hg19Gene(aDriverGene.gene());
        }
    }

    public boolean isValid(DriverGene driverGene) {
        return isValidPredicate.test(driverGene);
    }

    public String dndsGeneName(DriverGene driverGene) {
        return dndsGeneFunction.apply(driverGene);
    }

}
