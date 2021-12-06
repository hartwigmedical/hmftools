package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.gripss.common.SvData.hasLength;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SHORT_RESCUE_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.DEDUP;
import static com.hartwig.hmftools.gripss.filters.FilterType.PON;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.gripss.FilterCache;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.filters.FilterConstants;

public final class LinkRescue
{
    public static Set<Breakend> findRescuedBreakends(final LinkStore linkStore, final FilterCache filterCache, boolean rescueShortSVs)
    {
        Set<Breakend> rescuedBreakends = Sets.newHashSet();

        for(Breakend breakend : linkStore.getBreakendLinksMap().keySet())
        {
            if(rescuedBreakends.contains(breakend))
                continue;

            if(!filterCache.hasFilters(breakend))
                continue;

            if(!isRescueCandidate(breakend, filterCache, rescueShortSVs))
                continue;

            // this is a candidate for link-based rescue

            Set<Breakend> linkedBreakends = findLinkVariants(breakend, linkStore);

            if(linkedBreakends.isEmpty())
                continue;

            boolean hasValidPassing = linkedBreakends.stream()
                    .anyMatch(x -> !filterCache.hasFilters(x) && isRescueCandidate(x, filterCache, rescueShortSVs));

            if(hasValidPassing)
            {
                for(Breakend linkedBreakend : linkedBreakends)
                {
                    if(!isRescueCandidate(linkedBreakend, filterCache, rescueShortSVs))
                        continue;

                    rescuedBreakends.add(linkedBreakend);

                    Breakend otherLinkedBreakend = linkedBreakend.otherBreakend();

                    if(otherLinkedBreakend != null)
                        rescuedBreakends.add(otherLinkedBreakend);
                }
            }
        }

        return rescuedBreakends;
    }

    private static boolean isRescueCandidate(final Breakend breakend, final FilterCache filterCache, boolean rescueShortSVs)
    {
        if(filterCache.hasFilter(breakend.sv(), DEDUP))
            return false;

        if(!rescueShortSVs && hasLength(breakend.type()) && breakend.sv().length() < SHORT_RESCUE_LENGTH)
            return false;

        return true;
    }

    private static Set<Breakend> findLinkVariants(final Breakend breakend, final LinkStore linkStore)
    {
        Set<Breakend> linkedBreakends = Sets.newHashSet();
        findLinkVariants(breakend, linkStore, linkedBreakends);
        return linkedBreakends;
    }

    private static void findLinkVariants(final Breakend breakend, final LinkStore linkStore, final Set<Breakend> linkedBreakends)
    {
        if(linkedBreakends.contains(breakend))
            return;

        linkedBreakends.add(breakend);

        // call recursively on all other breakends that this one is linked to, and their SV's other breakend
        List<Link> links = linkStore.getBreakendLinks(breakend);

        if(links == null)
            return;

        for(Link link : links)
        {
            Breakend linkedBreakend = link.otherBreakend(breakend);
            findLinkVariants(linkedBreakend, linkStore, linkedBreakends);

            if(linkedBreakend.otherBreakend() != null)
                findLinkVariants(linkedBreakend.otherBreakend(), linkStore, linkedBreakends);
        }
    }

    public static Set<Breakend> findRescuedDsbLineInsertions(final LinkStore linkStore, final FilterCache filterCache, double minQual)
    {
        Set<Breakend> rescuedBreakends = Sets.newHashSet();

        for(Breakend breakend : linkStore.getBreakendLinksMap().keySet())
        {
            if(rescuedBreakends.contains(breakend))
                continue;

            if(!filterCache.hasFilters(breakend))
                continue;

            if(!isLineRescueCandidate(breakend, filterCache))
                continue;

            List<Link> links = linkStore.getBreakendLinks(breakend);

            for(Link link : links)
            {
                Breakend otherBreakend = link.otherBreakend(breakend);

                if(!breakend.IsLineInsertion && !otherBreakend.IsLineInsertion)
                    continue;

                double combinedQual = breakend.Qual + otherBreakend.Qual;

                if(combinedQual < minQual)
                    continue;

                if(!isLineRescueCandidate(otherBreakend, filterCache))
                    continue;

                rescuedBreakends.add(breakend);

                if(breakend.otherBreakend() != null)
                    rescuedBreakends.add(breakend.otherBreakend());

                rescuedBreakends.add(otherBreakend);

                if(otherBreakend.otherBreakend() != null)
                    rescuedBreakends.add(otherBreakend.otherBreakend());
            }
        }

        return rescuedBreakends;
    }

    private static boolean isLineRescueCandidate(final Breakend breakend, final FilterCache filterCache)
    {
        if(filterCache.hasFilter(breakend.sv(), PON))
            return false;

        if(filterCache.hasFilter(breakend.sv(), DEDUP))
            return false;

        if(hasLength(breakend.type()) && breakend.sv().length() < FilterConstants.SHORT_RESCUE_LENGTH)
            return false;

        return true;
    }

        /*
        data class LinkRescue(val rescues: Set<String>) {
        companion object {

            fun rescueDsb(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore): LinkRescue {
                return rescue(links, softFilterStore, variantStore, false)
            }

            fun rescueAssembly(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore): LinkRescue {
                return rescue(links, softFilterStore, variantStore, true)
            }

            fun rescueTransitive(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore): LinkRescue {
                return rescue(links, softFilterStore, variantStore, true)
            }

        private fun rescue(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore, rescueShort: Boolean): LinkRescue {
            val rescues = mutableSetOf<String>()
            val isValid = { x: StructuralVariantContext -> (rescueShort || !x.isTooShortToRescue) && !softFilterStore.containsDuplicateFilter(x.vcfId, x.mateId) }
            val isRescueCandidate = { x: StructuralVariantContext -> isValid(x) && softFilterStore.isFiltered(x.vcfId) }

            for (vcfId in links.linkedVariants()) {
                if (rescues.contains(vcfId)) {
                    continue
                }

                val variant = variantStore.select(vcfId)
                if (isRescueCandidate(variant)) {
                    val allLinkedVariants = allLinkedVariants(variant.vcfId, links, variantStore).map { variantStore.select(it) }
                    val anyPassing = allLinkedVariants.any { x -> isValid(x) && softFilterStore.isPassing(x.vcfId) }
                    if (anyPassing) {
                        for (linkedVariant in allLinkedVariants ){
                            if (isRescueCandidate(linkedVariant)) {
                                rescues.add(linkedVariant.vcfId)
                                linkedVariant.mateId?.let { rescues.add(it) }
                            }
                        }
                    }
                }
            }

            return LinkRescue(rescues)
        }

        private fun allLinkedVariants(vcfId: String, links: LinkStore, variantStore: VariantStore): Set<String> {
        val result = mutableSetOf<String>()
        allLinkedVariantsInner(vcfId, links, variantStore, result)
        return result
    }

        private fun allLinkedVariantsInner(vcfId: String, links: LinkStore, variantStore: VariantStore, result: MutableSet<String>) {
            if (result.contains(vcfId)) {
                return
            }

            result.add(vcfId)
            val linkedVariants = links.linkedVariants(vcfId).map { variantStore.select(it.otherVcfId) }
            for (linkedVariant in linkedVariants) {
                allLinkedVariantsInner(linkedVariant.vcfId, links, variantStore, result)
                linkedVariant.mateId?.let { allLinkedVariantsInner(it, links, variantStore, result) }
            }
        }

        fun rescueDsbMobileElementInsertion(config: GripssFilterConfig, links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore): LinkRescue {
            val rescues = mutableSetOf<String>()
            val isValid = { x: StructuralVariantContext -> !softFilterStore.containsPONFilter(x.vcfId, x.mateId) && !softFilterStore.containsDuplicateFilter(x.vcfId, x.mateId) && !x.isTooShortToRescue }
            val isRescueCandidate = { x: StructuralVariantContext -> softFilterStore.isFiltered(x.vcfId) && isValid(x) }

            for (vcfId in links.linkedVariants()) {
                if (rescues.contains(vcfId)) {
                    continue
                }

                for (link in links.linkedVariants(vcfId)) {
                    val variant = variantStore.select(vcfId)
                    if (isRescueCandidate(variant)) {
                        val other = variantStore.select(link.otherVcfId)
                        val combinedQual = variant.tumorQual + other.tumorQual
                        if (combinedQual >= config.minQualRescueMobileElementInsertion && isValid(other) && (variant.isMobileElementInsertion || other.isMobileElementInsertion)) {
                            rescues.add(variant.vcfId)
                            variant.mateId?.let { rescues.add(it) }
                            rescues.add(other.vcfId)
                            other.mateId?.let { rescues.add(it) }
                        }
                    }
                }
            }

            return LinkRescue(rescues)
        }
     */

}

