package com.hartwig.hmftools.sage.phase;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SimpleVariant;

public class LpsReadCounts
{
    public final int Id;
    public final List<SimpleVariant> Variants;
    public List<LpsReadCounts> SubsetVariantCounts;

    public int Depth;
    public int AltSupport;

    public LpsReadCounts(final int id, final List<SimpleVariant> variants)
    {
        Id = id;
        Variants = variants;
        Depth = 0;
        AltSupport = 0;
        SubsetVariantCounts = null;
    }

    public void registerCount(final List<SimpleVariant> variants, boolean hasAltSupport)
    {
        LpsReadCounts counts = this;

        if(variants.size() < Variants.size())
        {
            if(SubsetVariantCounts == null)
            {
                SubsetVariantCounts = Lists.newArrayList();
            }

            counts = SubsetVariantCounts.stream().filter(x -> x.variantsMatch(variants)).findFirst().orElse(null);

            if(counts == null)
            {
                counts = new LpsReadCounts(SubsetVariantCounts.size(), variants);
                SubsetVariantCounts.add(counts);
            }
        }

        ++counts.Depth;

        if(hasAltSupport)
        {
            ++counts.AltSupport;
        }
    }

    public boolean variantsMatch(final List<SimpleVariant> variants)
    {
        if(Variants.size() != variants.size())
            return false;

        return Variants.stream().allMatch(x -> variants.contains(x));
    }

    public String toString()
    {
        return format("%d: variants(%d) depth(%d) alt(%d) subGroups(%d)",
                Id, Variants.size(), Depth, AltSupport, SubsetVariantCounts != null ? SubsetVariantCounts.size() : 0);
    }
}
