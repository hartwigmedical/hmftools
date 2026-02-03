package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

public class VariantCache
{
    private final List<VariantPonData> mVariants;
    private int mLastVariantIndex; // to speed up searching

    public VariantCache()
    {
        mVariants = Lists.newArrayList();
        mLastVariantIndex = 0;
    }

    public VariantPonData   getOrCreateVariant(final String chromosome, final int position, final String ref, final String alt)
    {
        // start from the last inserted index since each VCF is ordered

        if(!mVariants.isEmpty() && mLastVariantIndex > 0 && position < mVariants.get(mLastVariantIndex).Position)
        {
            PV_LOGGER.error("variant cache invalid: size({}) lastIndex({} pos={}) new position({})",
                    mVariants.size(), mLastVariantIndex, mVariants.get(mLastVariantIndex).Position, position);

            mLastVariantIndex = 0; // reset so can recover
        }

        int index = mLastVariantIndex;
        while(index < mVariants.size())
        {
            VariantPonData variant = mVariants.get(index);

            if(position > variant.Position)
            {
                ++index;
                continue;
            }

            if(position < variant.Position)
                break;

            if(variant.Ref.equals(ref) && variant.Alt.equals(alt))
                return variant;

            ++index;
        }

        VariantPonData newVariant = new VariantPonData(chromosome, position, ref, alt);
        mVariants.add(index, newVariant);

        mLastVariantIndex = index;

        return newVariant;
    }

    public void resetSearch() { mLastVariantIndex = 0; }

    public void replaceVariants(final List<VariantPonData> variants)
    {
        mVariants.clear();
        mVariants.addAll(variants);
    }

    public int variantCount() { return mVariants.size(); }
    public List<VariantPonData> variants() { return mVariants; }

    @VisibleForTesting
    public int lastVariantIndex() { return mLastVariantIndex; }
}
