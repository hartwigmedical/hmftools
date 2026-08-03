package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;

import java.util.List;

import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class GeneCopyNumberData extends ComparableItem
{
    public final GeneCopyNumber CopyNumber;

    public GeneCopyNumberData(final GeneCopyNumber copyNumber, final List<FieldInfo> fields)
    {
        CopyNumber = copyNumber;

        addDoubleValue(GeneCopyNumberComparer.Fields.MinCopyNumber.toString(), copyNumber.MinCopyNumber, fields);
        addDoubleValue(GeneCopyNumberComparer.Fields.MaxCopyNumber.toString(), copyNumber.MaxCopyNumber, fields);

        // TODO: display only
        addIntValue(GeneCopyNumberComparer.Fields.MinRegionStart.toString(), copyNumber.MinRegionStart, fields);
        addIntValue(GeneCopyNumberComparer.Fields.MinRegionEnd.toString(), copyNumber.MinRegionEnd, fields);
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
