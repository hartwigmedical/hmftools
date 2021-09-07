package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.genome.genepanel.GeneNameMapping;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.jetbrains.annotations.NotNull;

class DndsGeneName
{
    private final GeneNameMapping mGeneNameMapping;
    private final Predicate<DriverGene> mAlignmentMapPredicate;
    private final Function<DriverGene, String> mDndsGeneFunction;
    private final Set<String> mDndsGenes;

    public DndsGeneName(@NotNull final RefGenomeVersion refGenomeVersion, @NotNull final Set<String> dndsGenes)
    {
        mDndsGenes = dndsGenes;
        mGeneNameMapping = GeneNameMapping.loadFromEmbeddedResource();

        if(refGenomeVersion == RefGenomeVersion.V37)
        {
            mAlignmentMapPredicate = driverGene -> mGeneNameMapping.isValidV37Gene(driverGene.gene());
            mDndsGeneFunction = DriverGene::gene;
        }
        else
        {
            mAlignmentMapPredicate = driverGene -> mGeneNameMapping.isValidV38Gene(driverGene.gene());
            mDndsGeneFunction = driverGene -> mGeneNameMapping.v37Gene(driverGene.gene());
        }
    }

    public boolean isValid(@NotNull DriverGene driverGene)
    {
        if(!mAlignmentMapPredicate.test(driverGene))
        {
            return false;
        }

        String dndsGeneName = dndsGeneName(driverGene);
        return mDndsGenes.contains(dndsGeneName);
    }

    @NotNull
    public String dndsGeneName(@NotNull DriverGene driverGene)
    {
        return mDndsGeneFunction.apply(driverGene);
    }
}
