package com.hartwig.hmftools.lilac.coverage;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class HlaComplex
{
    public final List<HlaAllele> Alleles;

    public HlaComplex(final List<HlaAllele> alleles)
    {
        Alleles = alleles;
    }

    public String toString()
    {
        return String.format("alleles(%s)", HlaAllele.toString(Alleles));
    }

    public static Set<HlaComplex> findDuplicates(final List<HlaComplex> complexes)
    {
        Set<HlaComplex> duplicates = Sets.newHashSet();

        for(int i = 0; i < complexes.size() - 1; ++i)
        {
            HlaComplex comp1 = complexes.get(i);

            for(int j = i + 1; j < complexes.size(); ++j)
            {
                if(comp1.matches(complexes.get(j)))
                {
                    duplicates.add(complexes.get(i));
                    break;
                }
            }
        }

        return duplicates;
    }

    public boolean matches(final HlaComplex other)
    {
        if(Alleles.size() != other.Alleles.size())
            return false;

        for(int i = 0; i < Alleles.size(); ++i)
        {
            if(!Alleles.get(i).equals(other.Alleles.get(i)))
                return false;
        }

        return true;
    }

}
