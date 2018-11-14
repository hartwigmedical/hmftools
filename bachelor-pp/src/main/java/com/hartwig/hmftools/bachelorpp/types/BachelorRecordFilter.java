package com.hartwig.hmftools.bachelorpp.types;

import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;

import javafx.geometry.Pos;

public class BachelorRecordFilter
{
    public final String Chromosome;
    public final long Position;
    public final String Id;
    public final String Ref;
    public final String Alt;
    public final String Qual;
    public final String Filter;
    public final String Info;

    public BachelorRecordFilter(final String chromosome, final long position, final String id, final String ref, final String alt,
            final String qual, final String filter, final String info)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Id = id;
        Qual = qual;
        Info = info;
        Filter = filter;
    }

    public boolean matches(final BachelorGermlineVariant var)
    {
        if(!var.chromosome().equals(Chromosome) || var.position() != Position
        || !var.ref().equals(Ref) || !var.alts().equals(Alt))
        {
            return false;
        }

        if(var.effects().equals(STOP_GAINED) && Info.contains("nonsense"))
            return true;

        if(!Info.contains(var.effects()))
            return false;

        return true;
    }

}
