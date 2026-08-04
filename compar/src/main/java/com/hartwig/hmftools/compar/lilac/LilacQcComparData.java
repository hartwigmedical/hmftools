package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.compar.common.CategoryType.LILAC_QC;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class LilacQcComparData extends ComparableItem
{
    public final com.hartwig.hmftools.common.hla.LilacQcData QcData;
    public final List<LilacAllele> Alleles;

    public LilacQcComparData(final LilacQcData qcData, final List<LilacAllele> alleles, final List<FieldInfo> fields)
    {
        QcData = qcData;
        Alleles = alleles;

        addStringValue(LilacQcComparer.Fields.Status.toString(), QcData.status(), fields);
        addIntValue(LilacQcComparer.Fields.TotalFragments.toString(), QcData.totalFragments(), fields);
        addIntValue(LilacQcComparer.Fields.FittedFragments.toString(), QcData.fittedFragments(), fields);
        addIntValue(LilacQcComparer.Fields.DiscardedAlignmentFragments.toString(), QcData.discardedAlignmentFragments(), fields);
        addIntValue(LilacQcComparer.Fields.DiscardedIndels.toString(), QcData.discardedIndels(), fields);
        addStringValue(LilacQcComparer.Fields.HlaYAllele.toString(), QcData.hlaYAllele(), fields);

        String allelesStr = Alleles.stream().map(x -> x.allele()).collect(Collectors.joining(ITEM_DELIM));
        addStringValue(LilacQcComparer.Fields.Alleles.toString(), allelesStr, fields);
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
        final LilacQcComparData otherData = (LilacQcComparData)other;

        if(!QcData.genes().equals(otherData.QcData.genes()))
            return false;

        return true;
    }
}
