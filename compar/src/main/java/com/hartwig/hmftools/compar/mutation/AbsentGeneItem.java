package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CategoryType.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.common.field.FieldInfo.findField;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_GENE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.FieldValue;
import com.hartwig.hmftools.compar.common.field.StringFieldValue;

public class AbsentGeneItem extends ComparableItem
{
    private final String mGeneName;
    private final FieldValue mGeneValue;

    public AbsentGeneItem(final String geneName, final List<FieldInfo> fields)
    {
        mGeneName = geneName;
        mGeneValue = new StringFieldValue(findField(FLD_GENE, fields), geneName);
    }

    @Override
    public CategoryType category() { return SOMATIC_VARIANT; }

    @Override
    public String key() { return mGeneName; }

    @Override
    public String geneName() { return mGeneName; }

    @Override
    public Map<String, FieldValue> fieldValues() { return Maps.newHashMap(Map.of(FLD_GENE, mGeneValue)); }

    @Override
    public boolean matches(final ComparableItem other) { return mGeneName.equals(other.geneName()); }
}
