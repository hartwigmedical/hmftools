package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.gripss.links.TransitiveLink.TRANS_LINK_PREFIX;

import java.util.ArrayDeque;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.SvDataCache;

public class TransitiveLinkFinder
{
    private final SvDataCache mSvDataCache;
    private final LinkStore mAssemblyLinkStore;

    // TODO: make public and move to Constants
    public static final int MAX_ASSEMBLY_JUMPS = 5;
    public static final int MAX_TRANSITIVE_JUMPS = 2;
    private static final int MAX_VARIANTS = 500000;
    private static final int MAX_ALTERNATIVES = 25;
    private static final int MAX_ALTERNATIVES_SEEK_DISTANCE = 1000;
    private static final int MAX_ALTERNATIVES_ADDITIONAL_DISTANCE = MAX_ALTERNATIVES_SEEK_DISTANCE;

    private static final int MIN_TRANSITIVE_DISTANCE = 30;
    private static final int MAX_TRANSITIVE_SEEK_DISTANCE = 2000;
    private static final int MAX_TRANSITIVE_ADDITIONAL_DISTANCE = 1000;

    public TransitiveLinkFinder(final SvDataCache svDataCache, final LinkStore assemblyLinkStore)
    {
        mSvDataCache = svDataCache;
        mAssemblyLinkStore = assemblyLinkStore;
    }

    public List<Link> findTransitiveLinks(final Breakend breakend)
    {
        List<Link> links = Lists.newArrayList();

        if(mSvDataCache.getSvList().size() > MAX_VARIANTS)
            return links;

        if(breakend.isSgl())
            return links;

        Breakend target = breakend.otherBreakend();

        List<Breakend> alternatives = selectAlternatives(breakend);

        if(alternatives.isEmpty() || alternatives.size() > MAX_ALTERNATIVES)
            return links;

        // logger.debug("Examining ${alternativeStarts.size} alternative(s) to variant $target")

        String tranksLinkPrefix = String.format("%s_%s_",TRANS_LINK_PREFIX, breakend.VcfId);

        ArrayDeque assemblyTransLinks = new ArrayDeque<TransitiveLink>();
        ArrayDeque transLinks = new ArrayDeque<TransitiveLink>();
        ArrayDeque matchedTransLinks = new ArrayDeque<TransitiveLink>();

        for(Breakend alternative : alternatives)
        {
            Breakend alternativeOtherBreakend = alternative.otherBreakend();

            Link altLink = Link.from(alternative.sv());

            TransitiveLink transLink = new TransitiveLink(tranksLinkPrefix, alternative, alternativeOtherBreakend, Lists.newArrayList(altLink));

            assemblyTransLinks.add(transLink);
        }

        List<Link> assemblyLinks = findLinks(target, assemblyTransLinks, transLinks, matchedTransLinks);

        if(!assemblyLinks.isEmpty())
            return assemblyLinks;

        return links;
    }

    private List<Link> findLinks(
            final Breakend target, final ArrayDeque assemblyTransLinks, final ArrayDeque transLinks, final ArrayDeque matchedTransLinks)
    {
        if(transLinks.size() > 1)
        {
            // no result if we there is more than one transitive path (and no assembly path)
            return Lists.newArrayList();
        }

        if (assemblyTransLinks.isEmpty() && transLinks.isEmpty())
        {
            if (matchedTransLinks.size() == 1)
            {
                TransitiveLink node = (TransitiveLink)matchedTransLinks.pop();
                return node.links();
            }

            return Lists.newArrayList();
        }

        if(!assemblyTransLinks.isEmpty())
        {
            TransitiveLink node = (TransitiveLink)assemblyTransLinks.removeFirst();

            if(node.matchesTarget(target))
            {
                // return the first (breath-wise) completely assembled link
                return node.links();
            }

            List<TransitiveLink> newAssemblyTransLinks = createAssemblyNodes(node);
            newAssemblyTransLinks.forEach(x -> assemblyTransLinks.add(x));

            List<TransitiveLink> newTransLinks = createTransitiveNodes(node);
            newTransLinks.forEach(x -> transLinks.add(x));
        }
        else
        {
            TransitiveLink node = (TransitiveLink)transLinks.removeFirst();

            if(node.matchesTarget(target))
            {
                matchedTransLinks.add(node);
            }

            List<TransitiveLink> newAssemblyTransLinks = createAssemblyNodes(node);
            newAssemblyTransLinks.forEach(x -> transLinks.add(x));

            List<TransitiveLink> newTransLinks = createTransitiveNodes(node);
            newTransLinks.forEach(x -> transLinks.add(x));
        }

        return findLinks(target, assemblyTransLinks, transLinks, matchedTransLinks);
    }

    private List<Breakend> selectAlternatives(final Breakend breakend)
    {
        List<Breakend> closeBreakends = mSvDataCache.selectOthersNearby(
                breakend, MAX_ALTERNATIVES_ADDITIONAL_DISTANCE, MAX_ALTERNATIVES_SEEK_DISTANCE);

        List<Breakend> alternatives = Lists.newArrayList();

        for(Breakend closeBreakend : closeBreakends)
        {
            if(closeBreakend.imprecise() || closeBreakend.isSgl())
                continue;

            if(!TransitiveLink.isAlternative(breakend, closeBreakend))
                continue;

            // add by Qual descending
            int index = 0;
            while(index < alternatives.size())
            {
                if(closeBreakend.Qual > alternatives.get(index).Qual)
                    break;

                ++index;
            }

            alternatives.add(index, closeBreakend);
        }

        return alternatives;
    }

    public List<TransitiveLink> createAssemblyNodes(final TransitiveLink transLink)
    {
        List<TransitiveLink> transitiveLinks = Lists.newArrayList();

        // create assembled links first
        if(transLink.remainingAssemblyJumps() == 0)
            return transitiveLinks;

        // get assembly links for the end of the transitive link
        Breakend transBreakend = transLink.breakendEnd();
        List<Link> unfilteredAssemblyLinks = mAssemblyLinkStore.getBreakendLinks(transBreakend);

        List<Link> assemblyLinkedBreakends = Lists.newArrayList();

        if(unfilteredAssemblyLinks != null)
        {
            for(Link assemblyLink : unfilteredAssemblyLinks)
            {
                if(transLink.links().contains(assemblyLink))
                    continue;

                Breakend otherBreakend = assemblyLink.otherBreakend(transBreakend);

                // add in order of quality of the other breakend
                int index = 0;
                while(index < assemblyLinkedBreakends.size())
                {
                    Breakend otherLinkBreakend = assemblyLinkedBreakends.get(index).otherBreakend(transBreakend);

                    if(otherBreakend.Qual > otherLinkBreakend.Qual)
                        break;

                    ++index;
                }

                assemblyLinkedBreakends.add(index, assemblyLink);
            }
        }

        for(Link assemblyLink : assemblyLinkedBreakends)
        {
            Breakend pairedBreakend = assemblyLink.otherBreakend(transBreakend);

            if(pairedBreakend.isSgl() || pairedBreakend.imprecise())
                continue;

            Breakend pairedOtherBreakend = pairedBreakend.otherBreakend();

            List<Link> newLinks = Lists.newArrayList(transLink.links());
            newLinks.add(assemblyLink);
            newLinks.add(Link.from(pairedBreakend.sv()));

            transitiveLinks.add(new TransitiveLink(
                    transLink.prefix(), pairedBreakend, pairedOtherBreakend,
                    transLink.remainingAssemblyJumps() - 1, transLink.remainingTransitiveJumps(), newLinks));
        }

        /*
        val unfilteredAssemblyLinks = assemblyLinkStore.linkedVariants(end.vcfId)
        if (remainingAssemblyJumps > 0) {
            val assemblyLinkedVariants = unfilteredAssemblyLinks
                    .filter { x -> !links.contains(x) }
                    .map { x -> Pair(x, variantStore.select(x.otherVcfId)) }
                    .sortedByDescending { x -> x.second.tumorQual }

            for (linkPair in assemblyLinkedVariants) {
                val link = linkPair.first
                val linkedVariant = linkPair.second
                if (!linkedVariant.isSingle && !linkedVariant.imprecise) {
                    val linkedVariantMate = variantStore.select(linkedVariant.mateId!!)
                    result.add(Node(transLinkPrefix, maxTransitiveJumps, remainingAssemblyJumps - 1, remainingTransitiveJumps,
                    linkedVariant, linkedVariantMate, links + link + Link(linkedVariant)))
                }
            }
        }
        */

        return transitiveLinks;
    }

    public List<TransitiveLink> createTransitiveNodes(final TransitiveLink transLink)
    {
        List<TransitiveLink> transLinks = Lists.newArrayList();

        Breakend transBreakend = transLink.breakendEnd();
        List<Link> assemblyLinks = mAssemblyLinkStore.getBreakendLinks(transBreakend);

        // skip if this breakend has assembly links
        if(assemblyLinks != null && !assemblyLinks.isEmpty())
            return transLinks;

        if(transLink.remainingTransitiveJumps() == 0)
            return transLinks;

        // find nearby breakends which also have no entries in the assembly links store

        List<Breakend> closeBreakends = mSvDataCache.selectOthersNearby(
                transBreakend, MAX_TRANSITIVE_ADDITIONAL_DISTANCE, MAX_TRANSITIVE_SEEK_DISTANCE);

        for(Breakend otherBreakend : closeBreakends)
        {
            if(otherBreakend.imprecise() || otherBreakend.isSgl())
                continue;

            if(otherBreakend.Orientation == transBreakend.Orientation)
                continue;

            if(!areCandidateLink(transBreakend, otherBreakend))
                continue;

            Breakend pairedOtherBreakend = otherBreakend.otherBreakend();

            List<Link> newLinks = Lists.newArrayList(transLink.links());
            String linkPrefix = String.format("%s{%d-%d}", transLink.prefix(), MAX_TRANSITIVE_JUMPS, transLink.remainingTransitiveJumps());
            Link transitiveLink = Link.from(linkPrefix, transBreakend, otherBreakend);
            newLinks.add(transitiveLink);
            newLinks.add(Link.from(otherBreakend.sv()));

            transLinks.add(new TransitiveLink(
                    transLink.prefix(), otherBreakend, pairedOtherBreakend,
                    transLink.remainingAssemblyJumps(), transLink.remainingTransitiveJumps() - 1, newLinks));
        }

        return transLinks;
    }

    private static boolean areCandidateLink(final Breakend first, final Breakend second)
    {
        if(first.Orientation == second.Orientation)
            return false;

        if(first.posOrient())
        {
            // other breakend must be less by the min TI length
            // val leftFilter: SvFilter = { other -> other.maxStart <= variant.minStart - MIN_TRANSITIVE_DISTANCE }
            return second.maxPosition() <= first.minPosition() - MIN_TRANSITIVE_DISTANCE;
        }
        else
        {
            // val rightFilter: SvFilter = { other -> other.minStart >= variant.maxStart + MIN_TRANSITIVE_DISTANCE }
            return second.minPosition() >= first.maxPosition() + MIN_TRANSITIVE_DISTANCE;
        }
    }

            /*
        private fun VariantStore.selectTransitive(variant: StructuralVariantContext): Collection<StructuralVariantContext> {
            val leftFilter: SvFilter = { other -> other.maxStart <= variant.minStart - MIN_TRANSITIVE_DISTANCE }
            val rightFilter: SvFilter = { other -> other.minStart >= variant.maxStart + MIN_TRANSITIVE_DISTANCE }
            val directionFilter: SvFilter = if (variant.orientation == 1.toByte()) leftFilter else rightFilter
            val transitiveFilter: SvFilter = { other -> other.orientation != variant.orientation && !other.imprecise && !other.isSingle }

            return selectOthersNearby(variant, MAX_TRANSITIVE_ADDITIONAL_DISTANCE, MAX_TRANSITIVE_SEEK_DISTANCE) { x -> directionFilter(x) && transitiveFilter(x) }.sortByQualDesc()
        }
         */

        /*
        if (assemblyLinkStore.linkedVariants(end.vcfId).isNotEmpty()) {
            // Cannot have transitive links if there are assembly links
            return Collections.emptyList()
        }

        val result = mutableListOf<Node>()

        // unfilteredAssemblyLinks.isEmpty() &&
        if (remainingTransitiveJumps > 0) {
            val unlinkedFilter: SvFilter = { assemblyLinkStore.linkedVariants(it.vcfId).isEmpty() }
            val transitiveLinkedVariants = variantStore.selectTransitive(end).filter { x -> unlinkedFilter(x) }
            for (linkedVariant in transitiveLinkedVariants) {
                val linkedVariantMate = variantStore.select(linkedVariant.mateId!!)
                val link = Link("$transLinkPrefix${maxTransitiveJumps - remainingTransitiveJumps}", Pair(end, linkedVariant))
                result.add(Node(transLinkPrefix, maxTransitiveJumps, remainingAssemblyJumps, remainingTransitiveJumps - 1, linkedVariant, linkedVariantMate, links + link + Link(linkedVariant)))
            }
        }

        return result
            */

}
