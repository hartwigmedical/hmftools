package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.compar.common.CategoryType.LILAC_QC;

import java.util.List;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class LilacQcData implements ComparableItem
{
    public final com.hartwig.hmftools.common.hla.LilacQcData QcData;
    public final List<LilacAllele> Alleles;

    public LilacQcData(final com.hartwig.hmftools.common.hla.LilacQcData qcData, final List<LilacAllele> alleles)
    {
        QcData = qcData;
        Alleles = alleles;
    }

    @Override
    public CategoryType category() { return LILAC_QC; }

    @Override
    public String key()
    {
        return QcData.genes();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final LilacQcData otherData = (LilacQcData)other;

        if(!QcData.genes().equals(otherData.QcData.genes()))
            return false;

        return true;
    }
}
