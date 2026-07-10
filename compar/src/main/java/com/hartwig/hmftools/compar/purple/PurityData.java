package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class PurityData implements ComparableItem
{
    public final PurityContext Purity;

    public PurityData(final PurityContext purityContext)
    {
        Purity = purityContext;
    }

    public CategoryType category() {
        return PURITY;
    }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }
}
