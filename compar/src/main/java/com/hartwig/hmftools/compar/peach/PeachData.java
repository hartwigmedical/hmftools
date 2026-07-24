package com.hartwig.hmftools.compar.peach;

import static com.hartwig.hmftools.compar.common.CategoryType.PEACH;

import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public class PeachData implements ComparableItem
{
    public final PeachGenotype Genotype;

    public PeachData(final PeachGenotype genotype)
    {
        Genotype = genotype;
    }

    @Override
    public CategoryType category()
    {
        return PEACH;
    }

    @Override
    public String key()
    {
        return Genotype.gene() + " " + Genotype.allele();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final PeachData otherData = (PeachData) other;
        if(!Genotype.gene().equals(otherData.Genotype.gene()))
        {
            return false;
        }
        return Genotype.allele().equals(otherData.Genotype.allele());
    }
}
