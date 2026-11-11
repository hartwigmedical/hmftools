package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CategoryType.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_GENE;

import java.util.List;

import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class AbsentGeneItem extends ComparableItem
{
    private final String mGeneName;

    public AbsentGeneItem(final String geneName, final List<FieldInfo> fields)
    {
        mGeneName = geneName;
        addStringValue(FLD_GENE, geneName, fields);
    }

    @Override
    public CategoryType category() { return SOMATIC_VARIANT; }

    @Override
    public String key() { return mGeneName; }

    @Override
    public String geneName() { return mGeneName; }

    @Override
    public boolean matches(final ComparableItem other) { return mGeneName.equals(other.geneName()); }
}
