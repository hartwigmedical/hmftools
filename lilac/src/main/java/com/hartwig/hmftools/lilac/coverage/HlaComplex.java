package com.hartwig.hmftools.lilac.coverage;

import java.util.List;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class HlaComplex
{
    private final List<HlaAllele> Alleles;

    public HlaComplex(final List<HlaAllele> alleles)
    {
        Alleles = alleles;
    }

    public List<HlaAllele> getAlleles() { return Alleles; }

    public String toString()
    {
        return String.format("alleles(%s)", HlaAllele.toString(Alleles));
    }
}
