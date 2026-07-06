package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.compar.common.CategoryType.LILAC;

import java.util.List;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class LilacData implements ComparableItem
{
    public final LilacQcData QcData;
    public final List<LilacAllele> Alleles;

    protected static final String FLD_ALLELES = "Alleles";
    
    public LilacData(final LilacQcData qcData, final List<LilacAllele> alleles)
    {
        QcData = qcData;
        Alleles = alleles;
    }

    @Override
    public CategoryType category() { return LILAC; }

    @Override
    public String key()
    {
        return QcData.genes();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final LilacData otherData = (LilacData)other;

        if(!QcData.genes().equals(otherData.QcData.genes()))
            return false;

        return true;
    }
}
