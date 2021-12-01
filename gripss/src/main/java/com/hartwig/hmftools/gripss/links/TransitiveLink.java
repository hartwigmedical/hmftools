package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.links.TransitiveLinkFinder.MAX_ASSEMBLY_JUMPS;
import static com.hartwig.hmftools.gripss.links.TransitiveLinkFinder.MAX_TRANSITIVE_JUMPS;

import java.util.List;

import com.hartwig.hmftools.gripss.common.Breakend;

public class TransitiveLink // previously 'Node'
{
    private final String mPrefix;

    private int mRemainingAssemblyJumps;
    private int mRemainingTransitiveJumps;

    private final Breakend[] mBreakends;
    private final List<Link> mLinks;

    private int mMinDistance;
    private int mMaxDistance;

    public static final String TRANS_LINK_PREFIX = "trs";

    public TransitiveLink(final String prefix, final Breakend start, final Breakend end, final List<Link> links)
    {
        this(prefix, start, end, MAX_ASSEMBLY_JUMPS, MAX_TRANSITIVE_JUMPS, links);
    }

    public TransitiveLink(
            final String prefix, final Breakend start, final Breakend end, int assemblyJumps, int transitiveJumps, final List<Link> links)
    {
        mPrefix = prefix;
        mRemainingAssemblyJumps = assemblyJumps;
        mRemainingTransitiveJumps = transitiveJumps;

        mBreakends = new Breakend[] { start, end };
        mLinks = links;

        mMinDistance = 0;
        mMaxDistance = 0;

        for (Link link : links)
        {
            mMinDistance += link.minDistance();
            mMaxDistance += link.maxDistance();
        }
    }

    public String prefix() { return mPrefix; }
    public Breakend breakendStart() { return mBreakends[SE_START]; }
    public Breakend breakendEnd() { return mBreakends[SE_END]; }

    public List<Link> links() { return mLinks; }

    public int remainingAssemblyJumps() { return mRemainingAssemblyJumps; }
    public void decrementAssemblyJumps() { if(mRemainingAssemblyJumps > 0) --mRemainingAssemblyJumps; }

    public int remainingTransitiveJumps() { return mRemainingTransitiveJumps; }
    public void decrementTransitiveJumps() { if(mRemainingTransitiveJumps > 0) --mRemainingTransitiveJumps; }

    public boolean matchesTarget(final Breakend targetEnd)
    {
        if(!isAlternative(targetEnd, breakendEnd()))
            return false;

        if(!targetEnd.imprecise())
        {
            int targetDistance = targetEnd.insertSequenceLength() + targetEnd.sv().duplicationLength();

            if (targetDistance < mMinDistance || targetDistance > mMaxDistance)
                return false;
        }

        return true;
    }

    public String toString()
    {
        return String.format("%s breaks(%s - %s) links(%d)", mPrefix, breakendStart(),  breakendEnd(), mLinks.size());
    }

    /*
    fun transitiveNodes(assemblyLinkStore: LinkStore, variantStore: VariantStore): List<Node> {
    }

    private fun VariantStore.selectTransitive(variant: StructuralVariantContext): Collection<StructuralVariantContext> {
        val leftFilter: SvFilter = { other -> other.maxStart <= variant.minStart - MIN_TRANSITIVE_DISTANCE }
        val rightFilter: SvFilter = { other -> other.minStart >= variant.maxStart + MIN_TRANSITIVE_DISTANCE }
        val directionFilter: SvFilter = if (variant.orientation == 1.toByte()) leftFilter else rightFilter
        val transitiveFilter: SvFilter = { other -> other.orientation != variant.orientation && !other.imprecise && !other.isSingle }

        return selectOthersNearby(variant, MAX_TRANSITIVE_ADDITIONAL_DISTANCE, MAX_TRANSITIVE_SEEK_DISTANCE) { x -> directionFilter(x) && transitiveFilter(x) }.sortByQualDesc()
    }
     */

    public static boolean isAlternative(final Breakend target, final Breakend other)
    {
        return isAlternative(target, other, 1);
    }

    public static boolean isAlternative(final Breakend target, final Breakend other, int additionalAllowance)
    {
        if(target == other || target.sv() == other.sv())
            return false;

        if(target.Orientation != other.Orientation)
            return false;

        int targetMinStart = target.minPosition();
        int targetMaxStart = target.maxPosition();

        if(target.posOrient())
        {
            targetMinStart -= additionalAllowance;
            targetMaxStart += target.insertSequenceLength() + additionalAllowance;
        }
        else
        {
            targetMinStart -= target.insertSequenceLength() + additionalAllowance;
            targetMaxStart += additionalAllowance;
        }

        int otherMinStart = other.minPosition();
        int otherMaxStart = other.maxPosition();

        if(other.posOrient())
        {
            otherMaxStart += other.insertSequenceLength();
        }
        else
        {
            otherMinStart -= other.insertSequenceLength();
        }

        /*
        public static int insertSequenceAdditionalDistance(variant: StructuralVariantContext): Pair<Int, Int> {
        return if (variant.orientation == 1.toByte()) {
            Pair(0, variant.insertSequenceLength)
        } else {
            Pair(variant.insertSequenceLength, 0)
        }

        int targetInsDistance = insertSequenceAdditionalDistance(target)
        val minStart = target.minStart - targetInsDistance.first - additionalAllowance
        val maxStart = target.maxStart + targetInsDistance.second + additionalAllowance

        val otherInsDistance = insertSequenceAdditionalDistance(other)
        val otherMinStart = other.minStart - otherInsDistance.first
        val otherMaxStart = other.maxStart + otherInsDistance.second
        */

        return otherMinStart <= targetMaxStart && otherMaxStart >= targetMinStart;
    }
}
