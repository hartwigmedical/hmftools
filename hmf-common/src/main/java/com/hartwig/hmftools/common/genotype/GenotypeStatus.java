package com.hartwig.hmftools.common.genotype;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;

public enum GenotypeStatus
{
    HOM_REF("HOM"),
    HET("HET"),
    HOM_ALT("HOM"),
    UNKNOWN("UNKNOWN");

    @NotNull
    private final String simpleDisplay;

    GenotypeStatus(@NotNull final String simpleDisplay)
    {
        this.simpleDisplay = simpleDisplay;
    }

    @NotNull
    public String simplifiedDisplay()
    {
        return simpleDisplay;
    }

    @NotNull
    public static GenotypeStatus fromGenotype(@NotNull Genotype genotype)
    {
        switch(genotype.getType())
        {
            case HET:
                return HET;
            case HOM_REF:
                return HOM_REF;
            case HOM_VAR:
                return HOM_ALT;
        }

        return UNKNOWN;
    }
}
