package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_CODING_MUTATION;
import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_MUTATION;
import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_SV;
import static com.hartwig.hmftools.common.wisp.CategoryType.SUBCLONAL_MUTATION;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_PROBES;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.NONREPORTABLE_SV_COUNT;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.SUBCLONAL_COUNT;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.wisp.CategoryType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class VariantSelector
{
    private static final Logger LOGGER = LogManager.getLogger(VariantSelector.class);

    public static List<Variant> selectVariants(List<Variant> variants)
    {
        // TODO: should have this strategy:
        //  1. order variants as is done now
        //  2. filter variants by non-probe-related checks
        //  3. if decide variant can be selected, generate its probe.
        //  4. evaluate probe. if it's good, select the variant and add the probe. if not, move to next variant

        variants.sort(new VariantComparator());

        List<Variant> selectedVariants = Lists.newArrayList();

        ProximateLocations registeredLocations = new ProximateLocations();
        Map<String, Integer> geneDisruptions = Maps.newHashMap();
        int[] typeCounts = new int[CategoryType.values().length];

        variants.forEach(variant -> variant.setSelectionStatus(SelectionStatus.NOT_SET));
        selectVariants(variants, selectedVariants, registeredLocations, geneDisruptions, typeCounts, true);
        selectVariants(variants, selectedVariants, registeredLocations, geneDisruptions, typeCounts, false);

        StringJoiner sj = new StringJoiner(" ");
        for(CategoryType type : CategoryType.values())
        {
            sj.add(format("%s=%d", type, typeCounts[type.ordinal()]));
        }
        LOGGER.debug("selected variant type counts: {}", sj.toString());

        return selectedVariants;
    }

    private static void selectVariants(final List<Variant> variants, List<Variant> selectedVariants, ProximateLocations registeredLocations,
            Map<String, Integer> geneDisruptions, int[] typeCounts, boolean firstPass)
    {
        for(Variant variant : variants)
        {
            if(selectedVariants.size() >= SAMPLE_PROBES)
            {
                break;
            }

            if(variant.isSelected())
            {
                continue;
            }

            boolean canSelect;

            if(variant.checkFilters())
            {
                if(firstPass)
                {
                    if(!variant.checkAndRegisterGeneLocation(geneDisruptions))
                    {
                        variant.setSelectionStatus(SelectionStatus.GENE_LOCATIONS);
                        canSelect = false;
                    }
                    else if(exceedsMaxByType(variant.categoryType(), OTHER_SV, typeCounts, NONREPORTABLE_SV_COUNT)
                            || exceedsMaxByType(variant.categoryType(), SUBCLONAL_MUTATION, typeCounts, SUBCLONAL_COUNT))
                    {
                        variant.setSelectionStatus(SelectionStatus.EXCEEDS_COUNT);
                        canSelect = false;
                    }
                    else if(!variant.passNonReportableFilters(true))
                    {
                        // This will be checked again on the second pass.
                        variant.setSelectionStatus(SelectionStatus.FILTERED);
                        canSelect = false;
                    }
                    else
                    {
                        canSelect = true;
                    }
                }
                else
                {
                    // If the variant failed the strict filters before, give it a chance to pass the relaxed filters now.
                    // Otherwise, the variant failed other checks and it can't be selected.
                    if(variant.selectionStatus() == SelectionStatus.FILTERED)
                    {
                        canSelect = variant.passNonReportableFilters(false);
                    }
                    else
                    {
                        canSelect = false;
                    }
                }
            }
            else
            {
                canSelect = true;
            }

            if(canSelect)
            {
                if(variant.checkAndRegisterLocation(registeredLocations))
                {
                    selectVariant(variant, selectedVariants, typeCounts);
                }
                else
                {
                    variant.setSelectionStatus(SelectionStatus.PROXIMATE);
                }
            }
        }
    }

    private static void selectVariant(Variant variant, List<Variant> selectedVariants, int[] typeCounts)
    {
        selectedVariants.add(variant);
        variant.setSelectionStatus(SelectionStatus.SELECTED);
        ++typeCounts[variant.categoryType().ordinal()];
    }

    private static boolean exceedsMaxByType(final CategoryType variantCategory, final CategoryType categoryType, final int[] typeCounts,
            int maxCount)
    {
        if(variantCategory != categoryType)
        {
            return false;
        }
        if(maxCount < 0)
        {
            return true;
        }
        return typeCounts[categoryType.ordinal()] >= maxCount;
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

            if(first.categoryType() == OTHER_MUTATION || first.categoryType() == OTHER_CODING_MUTATION
                    || first.categoryType() == SUBCLONAL_MUTATION)
            {
                // to randomise selection of the lowest priority variants (and avoid multiple selections from highly amplified regions)
                // use the inverse of position as the final comparison
                int locationHash1 = ((SomaticMutation) first).locationHash();
                int locationHash2 = ((SomaticMutation) second).locationHash();

                if(locationHash1 != locationHash2)
                {
                    return locationHash1 < locationHash2 ? -1 : 1;
                }
            }
            else if(first.copyNumber() != second.copyNumber())
            {
                return first.copyNumber() > second.copyNumber() ? -1 : 1;
            }

            return 0;
        }
    }
}
