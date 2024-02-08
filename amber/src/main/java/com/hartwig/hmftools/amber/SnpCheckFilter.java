package com.hartwig.hmftools.amber;

import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class SnpCheckFilter implements Predicate<PositionEvidence>
{
    private final Map<String,List<AmberSite>> mSnpLoci;

    public SnpCheckFilter(final Multimap<Chromosome, AmberSite> snpLoci)
    {
        mSnpLoci = Maps.newHashMap();

        for(AmberSite amberSite : snpLoci.values())
        {
            if(!amberSite.snpCheck())
                continue;

            List<AmberSite> chrSites = mSnpLoci.get(amberSite.Chromosome);

            if(chrSites == null)
            {
                chrSites = Lists.newArrayList();
                mSnpLoci.put(amberSite.Chromosome, chrSites);
            }

            chrSites.add(amberSite);
        }
    }

    @Override
    public boolean test(final PositionEvidence baseDepth)
    {
        List<AmberSite> chrSites = mSnpLoci.get(baseDepth.Chromosome);

        return chrSites != null && chrSites.stream()
                .anyMatch(x -> x.matches(baseDepth.Chromosome, baseDepth.Position, baseDepth.ref(), baseDepth.alt()));
    }
}
