package com.hartwig.hmftools.compar.peach;

import static com.hartwig.hmftools.compar.common.CategoryType.PEACH;

import java.util.List;

import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class PeachData extends ComparableItem
{
    public final PeachGenotype Genotype;

    public PeachData(final PeachGenotype genotype, final List<FieldInfo> fields)
    {
        Genotype = genotype;

        addIntValue(PeachComparer.Fields.AlleleCount.toString(), genotype.alleleCount(), fields);
        addStringValue(PeachComparer.Fields.Function.toString(), genotype.function(), fields);
        addStringValue(PeachComparer.Fields.Drugs.toString(), genotype.linkedDrugs(), fields);
        addStringValue(PeachComparer.Fields.PrescriptionUrls.toString(), genotype.urlPrescriptionInfo(), fields);
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
