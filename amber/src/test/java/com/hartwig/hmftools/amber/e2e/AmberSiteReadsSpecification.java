package com.hartwig.hmftools.amber.e2e;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

class AmberSiteReadsSpecification
{
    private final HumanChromosome Chromosome;
    private final int Index;
    private final int RefReadsCount;
    private final int AltReadsCount;
    private final int OtherReadsCount;

    public AmberSiteReadsSpecification(String label, String countsString)
    {
        String[] parts = label.split("_");
        Chromosome = HumanChromosome.fromString(parts[0]);
        Index = Integer.parseInt(parts[1]);

        String[] counts = countsString.split(",");
        RefReadsCount = Integer.parseInt(counts[0]);
        AltReadsCount = Integer.parseInt(counts[1]);
        OtherReadsCount = Integer.parseInt(counts[2]);
    }

    public Chromosome chromosome()
    {
        return Chromosome;
    }

    public int altCount()
    {
        return AltReadsCount;
    }

    public int refCount()
    {
        return RefReadsCount;
    }

    public int depth()
    {
        return RefReadsCount + AltReadsCount + OtherReadsCount;
    }

    public List<AmberSiteRead> resolve(ListMultimap<Chromosome, AmberSite> sites)
    {
        AmberSite site = sites.get(Chromosome).get(Index);
        List<AmberSiteRead> result = new ArrayList<>();
        for(int i = 0; i < RefReadsCount; i++)
        {
            result.add(new AmberSiteRead(site, AmberSiteRead.BaseInstruction.REF));
        }
        for(int i = 0; i < AltReadsCount; i++)
        {
            result.add(new AmberSiteRead(site, AmberSiteRead.BaseInstruction.ALT));
        }
        for(int i = 0; i < OtherReadsCount; i++)
        {
            result.add(new AmberSiteRead(site, AmberSiteRead.BaseInstruction.OTHER));
        }
        return result;
    }
}
