package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_INDEL;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_MISSENSE;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_NFS;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_REF_TOTAL;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_SPLICE;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_SYNON;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_TUMOR_CN;
import static com.hartwig.hmftools.common.hla.LilacAllele.FLD_TUMOR_TOTAL;
import static com.hartwig.hmftools.compar.common.CategoryType.LILAC_ALLELE;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.CommonUtils.findDiffs;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class LilacAlleleData implements ComparableItem
{
    public final LilacAllele Allele;
    public final int Index;  // indicates entry in case same haplotype is called twice

    public LilacAlleleData(LilacAllele allele, int index)
    {
        Allele = allele;
        Index = index;
    }

    @Override
    public CategoryType category() { return LILAC_ALLELE; }

    @Override
    public String key()
    {
        return String.format("%s:%s (%d)", Allele.genes(), Allele.allele(), Index);
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final LilacAlleleData otherData = (LilacAlleleData)other;

        if(!Allele.genes().equals(otherData.Allele.genes()) || !Allele.allele().equals(otherData.Allele.allele())
                || Index != otherData.Index)
        {
            return false;
        }

        return true;
    }
}
