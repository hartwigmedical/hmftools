package com.hartwig.hmftools.ctdna;

import static java.lang.String.format;

import static com.hartwig.hmftools.ctdna.CategoryType.OTHER_SV;
import static com.hartwig.hmftools.ctdna.CategoryType.REPORTABLE_MUTATION;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class VariantSelection
{
    public static List<Variant> selectVariants(final List<Variant> variants, final PvConfig config)
    {
        Collections.sort(variants, new VariantComparator());

        if(config.WriteAll)
            return variants;

        List<Variant> selectedVariants = Lists.newArrayList();

        int index = 0;
        int[] typeCounts = new int[CategoryType.values().length];

        while(selectedVariants.size() < config.ProbeCount && index < variants.size())
        {
            Variant variant = variants.get(index);

            if(variant.categoryType().ordinal() <= REPORTABLE_MUTATION.ordinal())
            {
                addVariant(variant, selectedVariants, typeCounts);
            }
            else if(variant.categoryType() == OTHER_SV)
            {
                if(typeCounts[OTHER_SV.ordinal()] < config.NonReportableSvCount)
                {
                    addVariant(variant, selectedVariants, typeCounts);
                }
            }
            else
            {
                addVariant(variant, selectedVariants, typeCounts);
            }

            ++index;
        }

        StringJoiner sj = new StringJoiner(" ");
        for(CategoryType type : CategoryType.values())
        {
            sj.add(format("%s=%d", type, typeCounts[type.ordinal()]));
        }

        PV_LOGGER.info("select variant type counts: {}", sj.toString());

        return selectedVariants;
    }

    private static void addVariant(final Variant variant, final List<Variant> selectedVariants, final int[] typeCounts)
    {
        selectedVariants.add(variant);
        ++typeCounts[variant.categoryType().ordinal()];
    }

    private static class VariantComparator implements Comparator<Variant>
    {
        public int compare(final Variant first, final Variant second)
        {
            // category
            // reported
            // VCN / copy number
            if(first.categoryType() != second.categoryType())
            {
                return first.categoryType().ordinal() < second.categoryType().ordinal() ? -1 : 1;
            }

            if(first.reported() != second.reported())
            {
                return first.reported() ? -1 : 1;
            }

            if(first.copyNumber() != second.copyNumber())
                return first.copyNumber() > second.copyNumber() ? -1 : 1;

            return 0;
        }
    }

}
