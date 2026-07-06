package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_VARIANT;

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class GermlineVariantData implements ComparableItem
{
    public final SmallVariant Variant;
    public final Set<String> Filters;
    public final BasePosition mComparisonPosition;

    public static final String FILTER_DELIMITER = ";";

    public GermlineVariantData(final SmallVariant variant, final BasePosition comparisonPosition)
    {
        Variant = variant;
        Filters = Arrays.stream(variant.filter().split(FILTER_DELIMITER, -1)).collect(Collectors.toSet());
        mComparisonPosition = comparisonPosition;
    }

    @Override
    public CategoryType category() { return GERMLINE_VARIANT; }

    @Override
    public String key()
    {
        if(mComparisonPosition.Position != Variant.position())
        {
            return String.format("%s:%d %s>%s %s liftover(%s)",
                    Variant.chromosome(), Variant.position(), Variant.ref(), Variant.alt(), Variant.type(), mComparisonPosition);
        }
        else
        {
            return String.format("%s:%d %s>%s %s", Variant.chromosome(), Variant.position(), Variant.ref(), Variant.alt(), Variant.type());
        }
    }

    @Override
    public boolean reportable()
    {
        return Variant.reported();
    }

    @Override
    public String geneName() { return Variant.gene(); }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GermlineVariantData otherVar = (GermlineVariantData) other;

        if(!mComparisonPosition.Chromosome.equals(otherVar.Variant.chromosome()) || mComparisonPosition.Position != otherVar.Variant.position())
            return false;

        if(!Variant.ref().equals(otherVar.Variant.ref()) || !Variant.alt().equals(otherVar.Variant.alt()))
            return false;

        if(Variant.type() != otherVar.Variant.type())
            return false;

        return true;
    }
}
