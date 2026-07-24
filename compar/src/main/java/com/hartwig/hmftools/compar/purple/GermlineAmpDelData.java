package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.purple.ReportedStatus.REPORTED;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_AMP_DEL;

import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class GermlineAmpDelData implements ComparableItem
{
    public final GermlineAmpDel AmpDelData;
    public final String mComparisonChromosome;

    public GermlineAmpDelData(final GermlineAmpDel germlineAmpDel, final String comparisonChromosome)
    {
        AmpDelData = germlineAmpDel;
        mComparisonChromosome = comparisonChromosome;
    }

    public CategoryType category() {
        return GERMLINE_AMP_DEL;
    }

    @Override
    public String key()
    {
        return format("%s", AmpDelData.GeneName);
    }

    @Override
    public boolean reportable() {
        return AmpDelData.Reported == REPORTED;
    }

    @Override
    public String geneName() { return AmpDelData.GeneName; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GermlineAmpDelData otherAmpDel = (GermlineAmpDelData)other;
        return AmpDelData.GeneName.equals(otherAmpDel.AmpDelData.GeneName);
    }
}
