package com.hartwig.hmftools.compar.sigs;

import java.util.List;

import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class SigsData extends ComparableItem
{
    public final SignatureAllocation Allocation;

    public SigsData(final SignatureAllocation allocation, final List<FieldInfo> fields)
    {
        Allocation = allocation;

        addDoubleValue(SigsComparer.Fields.Percent.toString(), allocation.percent(), fields);
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.SIGS;
    }

    @Override
    public String key()
    {
        return Allocation.signature();
    }

    @Override
    public boolean reportable()
    {
        return false;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final SigsData otherData = (SigsData)other;

        return otherData.Allocation.signature().equals(Allocation.signature());
    }
}
