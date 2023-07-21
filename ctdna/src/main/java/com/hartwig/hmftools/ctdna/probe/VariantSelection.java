package com.hartwig.hmftools.ctdna.probe;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.probe.CategoryType.OTHER_CODING_MUTATION;
import static com.hartwig.hmftools.ctdna.probe.CategoryType.OTHER_MUTATION;
import static com.hartwig.hmftools.ctdna.probe.CategoryType.OTHER_SV;
import static com.hartwig.hmftools.ctdna.probe.CategoryType.REPORTABLE_MUTATION;
import static com.hartwig.hmftools.ctdna.probe.CategoryType.SUBCLONAL_MUTATION;
import static com.hartwig.hmftools.ctdna.probe.SelectionStatus.EXCEEDS_COUNT;
import static com.hartwig.hmftools.ctdna.probe.SelectionStatus.FILTERED;
import static com.hartwig.hmftools.ctdna.probe.SelectionStatus.GENE_LOCATIONS;
import static com.hartwig.hmftools.ctdna.probe.SelectionStatus.PROXIMATE;
import static com.hartwig.hmftools.ctdna.probe.SelectionStatus.SELECTED;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public final class VariantSelection
{
    public static List<Variant> selectVariants(final List<Variant> variants, final ProbeConfig config)
    {
        Collections.sort(variants, new VariantComparator());

        List<Variant> selectedVariants = Lists.newArrayList();

        ProximateLocations registeredLocations = new ProximateLocations();
        Map<String,Integer> geneDisruptions = Maps.newHashMap();

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
            else
            {
                if(!variant.checkAndRegisterGeneLocation(geneDisruptions))
                {
                    canAdd = false;
                    variant.setSelectionStatus(GENE_LOCATIONS);
                }
                else if(exceedsMaxByType(variant.categoryType(), OTHER_SV, typeCounts, config.NonReportableSvCount)
                || exceedsMaxByType(variant.categoryType(), SUBCLONAL_MUTATION, typeCounts, config.SubclonalCount))
                {
                    canAdd = false;
                    variant.setSelectionStatus(EXCEEDS_COUNT);
                }
                else if(!variant.passNonReportableFilters(config))
                {
                    variant.setSelectionStatus(FILTERED);
                    canAdd = false;
                }
                else
                {
                    canAdd = true;
                }
            }

            if(canAdd)
            {
                if(variant.checkAndRegisterLocation(registeredLocations))
                {
                    addVariant(variant, selectedVariants, typeCounts);
                }
                else
                {
                    variant.setSelectionStatus(PROXIMATE);
                }
            }

            if(!variant.isSelected() && config.WriteAll)
                selectedVariants.add(variant);

            int variantSequences = selectedVariants.stream().filter(x -> x.isSelected()).mapToInt(x -> x.sequenceCount()).sum();
            if(variantSequences >= config.ProbeCount)
                break;

            ++index;
        }

        StringJoiner sj = new StringJoiner(" ");
        for(CategoryType type : CategoryType.values())
        {
            sj.add(format("%s=%d", type, typeCounts[type.ordinal()]));
        }

        CT_LOGGER.debug("selected variant type counts: {}", sj.toString());

        return selectedVariants;
    }

    private static boolean exceedsMaxByType(
            final CategoryType variantCategory, final CategoryType categoryType, final int[] typeCounts, final int maxCount)
    {
        if(variantCategory != categoryType)
            return false;

        if(maxCount < 0)
            return true;

        if(maxCount == 0)
            return false;

        return typeCounts[categoryType.ordinal()] >= maxCount;
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
        variant.setSelectionStatus(SELECTED);
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

            if(first.categoryType() == OTHER_MUTATION || first.categoryType() == OTHER_CODING_MUTATION || first.categoryType() == SUBCLONAL_MUTATION)
            {
                // to randomise selection of the lowest priority variants (and avoid multiple selections from highly amplified regions)
                // use the inverse of position as the final comparison
                int locationHash1 = ((SomaticMutation)first).locationHash();
                int locationHash2 = ((SomaticMutation)second).locationHash();

                if(locationHash1 != locationHash2)
                    return locationHash1 < locationHash2 ? -1 : 1;
            }
            else
            {
                if(first.copyNumber() != second.copyNumber())
                    return first.copyNumber() > second.copyNumber() ? -1 : 1;
            }

            return 0;
        }
    }
}
