package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.genome.refgenome.GeneNameMapping;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.jetbrains.annotations.NotNull;

class DndsGeneName {

    private final GeneNameMapping geneNameMapping;
    private final Predicate<DriverGene> alignmentMapPredicate;
    private final Function<DriverGene, String> dndsGeneFunction;
    private final Set<String> dndsGenes;

    public DndsGeneName(@NotNull final RefGenomeVersion refGenomeVersion, @NotNull final Set<String> dndsGenes) {
        this.dndsGenes = dndsGenes;
        this.geneNameMapping = new GeneNameMapping();
        if (refGenomeVersion == RefGenomeVersion.V37) {
            alignmentMapPredicate = aDriverGene -> geneNameMapping.isValidV37Gene(aDriverGene.gene());
            dndsGeneFunction = DriverGene::gene;
        } else {
            alignmentMapPredicate = aDriverGene -> geneNameMapping.isValidV38Gene(aDriverGene.gene());
            dndsGeneFunction = aDriverGene -> geneNameMapping.v37Gene(aDriverGene.gene());
        }
    }

    public boolean isValid(@NotNull DriverGene driverGene) {
        if (!alignmentMapPredicate.test(driverGene)) {
            return false;
        }

        String dndsGeneName = dndsGeneName(driverGene);
        return dndsGenes.contains(dndsGeneName);
    }

    public String dndsGeneName(@NotNull DriverGene driverGene) {
        return dndsGeneFunction.apply(driverGene);
    }
}
