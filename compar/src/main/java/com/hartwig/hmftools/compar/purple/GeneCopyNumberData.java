package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;

import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class GeneCopyNumberData implements ComparableItem
{
    public final GeneCopyNumber CopyNumber;

    public GeneCopyNumberData(final GeneCopyNumber copyNumber)
    {
        CopyNumber = copyNumber;
    }

    public CategoryType category() {
        return GENE_COPY_NUMBER;
    }

    @Override
    public String key()
    {
        return format("%s", CopyNumber.geneName());
    }

    @Override
    public boolean reportable() {
        return false;
    }

    @Override
    public String geneName() { return CopyNumber.GeneName; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GeneCopyNumberData otherCn = (GeneCopyNumberData)other;

        return CopyNumber.geneName().equals(otherCn.CopyNumber.geneName());
    }
}
