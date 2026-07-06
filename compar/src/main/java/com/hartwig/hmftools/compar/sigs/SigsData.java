package com.hartwig.hmftools.compar.sigs;

import static com.hartwig.hmftools.common.sigs.SignatureAllocationFile.PERCENT_FLD;

import static org.apache.commons.lang3.StringUtils.capitalize;

import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public record SigsData(SignatureAllocation SignatureAllocation) implements ComparableItem
{
    public static String FLD_PERCENT = capitalize(PERCENT_FLD);

    @Override
    public CategoryType category()
    {
        return CategoryType.SIGS;
    }

    @Override
    public String key()
    {
        return SignatureAllocation.signature();
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

        return otherData.SignatureAllocation.signature().equals(SignatureAllocation.signature());
    }
}
