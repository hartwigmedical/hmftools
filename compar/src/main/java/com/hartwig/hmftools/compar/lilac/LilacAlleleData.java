package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.compar.common.CategoryType.LILAC_ALLELE;

import java.util.List;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class LilacAlleleData extends ComparableItem
{
    public final LilacAllele Allele;
    public final int Index;  // indicates entry in case same haplotype is called twice

    public LilacAlleleData(LilacAllele allele, int index, final List<FieldInfo> fields)
    {
        Allele = allele;
        Index = index;

        addIntValue(LilacAlleleComparer.Fields.RefTotal.toString(), Allele.refFragments(), fields);
        addIntValue(LilacAlleleComparer.Fields.TumorTotal.toString(), Allele.tumorFragments(), fields);
        addDoubleValue(LilacAlleleComparer.Fields.TumorCopyNumber.toString(), Allele.tumorCopyNumber(), fields);
        addDoubleValue(LilacAlleleComparer.Fields.SomaticMissense.toString(), Allele.somaticMissense(), fields);
        addDoubleValue(LilacAlleleComparer.Fields.SomaticNonsenseOrFrameshift.toString(), Allele.somaticNonsenseOrFrameshift(), fields);
        addDoubleValue(LilacAlleleComparer.Fields.SomaticSplice.toString(), Allele.somaticSplice(), fields);
        addDoubleValue(LilacAlleleComparer.Fields.SomaticInframeIndel.toString(), Allele.somaticInframeIndel(), fields);
        addDoubleValue(LilacAlleleComparer.Fields.SomaticSynonymous.toString(), Allele.somaticSynonymous(), fields);
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
