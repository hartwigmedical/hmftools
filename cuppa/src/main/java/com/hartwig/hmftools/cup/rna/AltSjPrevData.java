package com.hartwig.hmftools.cup.rna;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class AltSjPrevData
{
    public final String GeneId;
    public final String Type;
    public final BaseRegion Location;

    // for cancer type data
    public final Map<String,Double> CancerPrevalences;
    public double PrevalenceTotal;
    public double MinPrevalence;

    // for sample data
    public double Prevalence;
    public final double RawPrevalence;

    public AltSjPrevData(final String gene, final String type, final String chromosome, int posStart, int posEnd)
    {
        GeneId = gene;
        Type = type;
        Location = new BaseRegion(chromosome, posStart, posEnd);
        RawPrevalence = 0;
        Prevalence = 0;

        CancerPrevalences = null;
        PrevalenceTotal = 0;
        MinPrevalence = 0;
    }

    public AltSjPrevData(final String gene, final String type,
            final String chromosome, int posStart, int posEnd, double minPrevalence)
    {
        GeneId = gene;
        Type = type;
        Location = new BaseRegion(chromosome, posStart, posEnd);
        RawPrevalence = 0;
        Prevalence = 0;

        CancerPrevalences = Maps.newHashMap();
        PrevalenceTotal = 0;
        MinPrevalence = minPrevalence;
    }

    public boolean matches(final AltSjPrevData other)
    {
        return Location.matches(other.Location); // && Type.equals(other.Type);
    }

    public String toString() { return String.format("gene(%s) type(%s) location(%s) prev(%.4f adj=%.4f)",
            GeneId, Type, Location, RawPrevalence, Prevalence); }

}
