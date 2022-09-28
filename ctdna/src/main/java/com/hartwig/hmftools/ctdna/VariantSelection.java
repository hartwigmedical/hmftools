package com.hartwig.hmftools.ctdna;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.ctdna.CategoryType.OTHER_SV;
import static com.hartwig.hmftools.ctdna.CategoryType.REPORTABLE_MUTATION;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.VectorUtils;

public final class VariantSelection
{
    public static List<Variant> selectVariants(final List<Variant> variants, final PvConfig config)
    {
        Collections.sort(variants, new VariantComparator());

        if(config.WriteAll)
            return variants;

        List<Variant> selectedVariants = Lists.newArrayList();

        Map<String,List<Integer>> registeredLocations = Maps.newHashMap();

        int index = 0;
        int[] typeCounts = new int[CategoryType.values().length];

        while(index < variants.size())
        {
            Variant variant = variants.get(index);

            boolean canAdd = false;

            if(variant.categoryType().ordinal() <= REPORTABLE_MUTATION.ordinal())
            {
                canAdd = true;
            }
            else if(variant.categoryType() == OTHER_SV)
            {
                if(typeCounts[OTHER_SV.ordinal()] < config.NonReportableSvCount)
                {
                    canAdd = true;
                }
            }
            else
            {
                canAdd = true;
            }

            if(canAdd)
            {
                if(variant.checkAndRegisterLocation(registeredLocations))
                {
                    addVariant(variant, selectedVariants, typeCounts);
                }
            }

            int variantSequences = selectedVariants.stream().mapToInt(x -> x.sequenceCount()).sum();
            if(variantSequences >= config.ProbeCount)
                break;

            ++index;
        }

        StringJoiner sj = new StringJoiner(" ");
        for(CategoryType type : CategoryType.values())
        {
            sj.add(format("%s=%d", type, typeCounts[type.ordinal()]));
        }

        PV_LOGGER.info("selected variant type counts: {}", sj.toString());

        return selectedVariants;
    }

    public static final int NEAR_DISTANCE = 50;

    public static boolean isNearRegisteredLocation(
            final Map<String,List<Integer>> registeredLocations, final String chromosome, final int position)
    {
        List<Integer> positions = registeredLocations.get(chromosome);

        if(positions == null)
            return false;

        return positions.stream().anyMatch(x -> abs(x - position) <= NEAR_DISTANCE);
    }

    public static void addRegisteredLocation(final Map<String,List<Integer>> registeredLocations, final String chromosome, final int position)
    {
        List<Integer> positions = registeredLocations.get(chromosome);

        if(positions == null)
        {
            positions = Lists.newArrayList(position);
            registeredLocations.put(chromosome, positions);
            return;
        }

        int index = 0;
        while(index < positions.size())
        {
            if(positions.get(index) > position)
                ++index;

            break;
        }

        positions.add(index, position);
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
