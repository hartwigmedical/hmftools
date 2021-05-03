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
        this.geneNameMapping = GeneNameMapping.loadFromEmbeddedResource();
        if (refGenomeVersion == RefGenomeVersion.V37) {
            alignmentMapPredicate = driverGene -> geneNameMapping.isValidV37Gene(driverGene.gene());
            dndsGeneFunction = DriverGene::gene;
        } else {
            alignmentMapPredicate = driverGene -> geneNameMapping.isValidV38Gene(driverGene.gene());
            dndsGeneFunction = driverGene -> geneNameMapping.v37Gene(driverGene.gene());
        }
    }

    public boolean isValid(@NotNull DriverGene driverGene) {
        if (!alignmentMapPredicate.test(driverGene)) {
            return false;
        }

        String dndsGeneName = dndsGeneName(driverGene);
        return dndsGenes.contains(dndsGeneName);
    }

    @NotNull
    public String dndsGeneName(@NotNull DriverGene driverGene) {
        return dndsGeneFunction.apply(driverGene);
    }
}
