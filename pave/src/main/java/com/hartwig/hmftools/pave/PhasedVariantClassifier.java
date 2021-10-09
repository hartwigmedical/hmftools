package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

public class PhasedVariantClassifier
{
    private final List<PhasedVariants> mPhasedVariants;

    public PhasedVariantClassifier()
    {
        mPhasedVariants = Lists.newArrayList();
    }

    public void checkAddVariant(final VariantData variant)
    {
        if(!variant.hasLocalPhaseSet())
            return;

        PhasedVariants phasedVariants = mPhasedVariants.stream()
                .filter(x -> x.LocalPhaseId == variant.localPhaseSet()).findFirst().orElse(null);

        if(phasedVariants == null)
        {
            phasedVariants = new PhasedVariants(variant.localPhaseSet());
            mPhasedVariants.add(phasedVariants);
        }

        phasedVariants.Variants.add(variant);
    }

    public boolean hasCompleteVariants(int currentLocalPhaseSet)
    {
        return !mPhasedVariants.isEmpty() && mPhasedVariants.get(0).LocalPhaseId != currentLocalPhaseSet;
    }

    public List<PhasedVariants> popCompletePhasedVariants(int currentLocalPhaseSet)
    {
        List<PhasedVariants> completeVariants = Lists.newArrayList();

        while(!mPhasedVariants.isEmpty())
        {
            if(mPhasedVariants.get(0).LocalPhaseId != currentLocalPhaseSet)
            {
                completeVariants.add(mPhasedVariants.get(0));
                mPhasedVariants.remove(0);
            }
            else
            {
                break;
            }
        }

        return completeVariants;
    }

    public List<PhasedVariants> allPhasedVariants() { return mPhasedVariants; }
    public void clear() { mPhasedVariants.clear(); }

    public static void reclassifyPhasedVariants(final PhasedVariants variants)
    {
        if(variants.Variants.size() < 2)
            return;

        // PV_LOGGER.debug("processing {} phased variant variants with localPhaseId({})", variants.Variants.size(), variants.LocalPhaseId);

    }

}
