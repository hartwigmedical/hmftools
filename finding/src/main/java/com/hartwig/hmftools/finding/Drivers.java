package com.hartwig.hmftools.finding;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;

import org.jetbrains.annotations.Nullable;

final class Drivers
{
    private static final Set<PurpleDriverType> MUTATION_DRIVER_TYPES =
            Set.of(PurpleDriverType.MUTATION, PurpleDriverType.GERMLINE_MUTATION);

    public static List<PurpleDriver> nonCanonicalMutationEntries(List<PurpleDriver> drivers)
    {
        List<PurpleDriver> nonCanonicalVariantEntries = new java.util.ArrayList<>();
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
    public static PurpleDriver canonicalMutationEntryForGene(List<PurpleDriver> drivers, String geneToFind)
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
