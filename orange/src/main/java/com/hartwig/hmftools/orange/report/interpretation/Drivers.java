package com.hartwig.hmftools.orange.report.interpretation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Drivers
{
    private static final Set<PurpleDriverType> MUTATION_DRIVER_TYPES =
            Sets.newHashSet(PurpleDriverType.MUTATION, PurpleDriverType.GERMLINE_MUTATION);

    public static List<PurpleDriver> nonCanonicalMutationEntries(final List<PurpleDriver> drivers)
    {
        List<PurpleDriver> nonCanonicalVariantEntries = Lists.newArrayList();
        for(PurpleDriver driver : drivers)
        {
            if(MUTATION_DRIVER_TYPES.contains(driver.type()) && !driver.isCanonical())
            {
                nonCanonicalVariantEntries.add(driver);
            }
        }
        return nonCanonicalVariantEntries;
    }

    @Nullable
    public static PurpleDriver canonicalMutationEntryForGene(final List<PurpleDriver> drivers, final String geneToFind)
    {
        PurpleDriver highest = null;
        for(PurpleDriver driver : drivers)
        {
            if(MUTATION_DRIVER_TYPES.contains(driver.type()) && driver.gene().equals(geneToFind) && driver.isCanonical())
            {
                if(highest == null || driver.driverLikelihood() > highest.driverLikelihood())
                {
                    highest = driver;
                }
            }
        }

        return highest;
    }
}
