package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.Q_ARM;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.types.LinxConstants.SHORT_DB_LENGTH;
import static com.hartwig.hmftools.linx.types.LinxConstants.SHORT_TI_LENGTH;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ClusterMetrics
{
    public int OriginArms;
    public int FragmentArms;
    public int ConsistentArms;
    public int ComplexArms;

    public int DBCount;
    public int ShortDBCount;

    // sum of deletion bridges including regions interrupted by other SVs (eg unphased)
    public long TotalDBLength;

    // where (typically long) TIs cross unclustered deleted regions (typically simple DELs)
    // an indication that the cluster could have included these in shattering events
    public int TraversedDelLength;
    public int TraversedDelCount;
    public int TraversedShortDelCount;
    public int ImpliedTICount; // from traversed DELs

    public long ChainedLength; // sum of chain lengths, allowing for traversal of the same genomic region (ie from a DUP)

    public long TotalRange; // chained distance without double-count and less any isolated short TIs
    public long TraversedRange; // bases covered by the cluster using ploidy data
    public long TotalDeleted; // bases with traversed range but without cluster CN support

    // output from chaining routine - % of copy number segments had discernable A, B and cluster ploidy
    public double ValidAlleleJcnSegmentPerc;

    public int IndelCount; // count of indels on TIs
    public double IndelProbability; // of seeing X indels within the range of the cluster

    public ClusterMetrics()
    {
        OriginArms = 0;
        FragmentArms = 0;
        ConsistentArms = 0;
        ComplexArms = 0;
        DBCount = 0;
        ShortDBCount = 0;
        TotalDBLength = 0;
        TraversedDelLength = 0;
        TraversedDelCount = 0;
        TraversedShortDelCount = 0;
        ImpliedTICount = 0;
        ChainedLength = 0;
        TotalRange = 0;
        TraversedRange = 0;
        TotalDeleted = 0;
        ValidAlleleJcnSegmentPerc = 1;
        IndelCount = 0;
        IndelProbability = 1;
    }

    public void populate(final SvCluster cluster, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        if(cluster.getSvCount() == 0)
            return;

        final Map<String,Boolean> crossCentromereLinks = Maps.newHashMap();

        cluster.getLinkedPairs().stream()
                .filter(x -> x.firstBreakend().arm() != x.secondBreakend().arm())
                .forEach(x -> crossCentromereLinks.put(x.chromosome(), true));

        for(final Map.Entry<String,List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();
            final List<SvBreakend> chrBreakendList = chrBreakendMap.get(chromosome);

            int startIndex = findStartIndex(breakendList);
            int endIndex = findEndIndex(breakendList);

            if(startIndex < endIndex)
            {
                TotalRange += breakendList.get(endIndex).position() - breakendList.get(startIndex).position();
            }

            // add up all stand-alone DELs and deletion bridges, including those spanning unclustered SVs
            for(int i = startIndex; i <= endIndex - 1; ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvVarData var = breakend.getSV();
                final SvBreakend nextBreakend = breakendList.get(i+1);
                int breakendsDistance = nextBreakend.position() - breakend.position();

                boolean isDB = breakend.getDBLink() != null && breakend.getDBLink() == nextBreakend.getDBLink();

                boolean isSimpleDel = !isDB && var.type() == DEL
                        && breakend.orientation() == 1 && nextBreakend.getSV() == var;

                boolean isTraversedDB = !isDB && !isSimpleDel && breakend.orientation() == 1 && nextBreakend.orientation() == -1;

                if(isDB || isSimpleDel || isTraversedDB)
                {
                    ++DBCount;

                    int delLength = isDB ? max(breakend.getDBLink().length(), 0) : (isSimpleDel ? var.length() : breakendsDistance);

                    TotalDBLength += delLength;

                    if(delLength <= SHORT_DB_LENGTH)
                        ++ShortDBCount;
                }
                else if(nextBreakend.getChrPosIndex() > breakend.getChrPosIndex() + 2 && breakend.orientation() == -1 && nextBreakend.orientation() == 1)
                {
                    // a gap in the cluster breakends - check for DELs from other clusters
                    int impliedTICount = 0;

                    for(int j = breakend.getChrPosIndex() + 1; j < nextBreakend.getChrPosIndex() - 1; ++j)
                    {
                        final SvBreakend ncBreakend = chrBreakendList.get(j);
                        final SvVarData ncVar = ncBreakend.getSV();
                        if(ncVar.type() == DEL && chrBreakendList.get(j + 1).getSV() == ncVar)
                        {
                            TraversedDelLength += ncVar.length();
                            ++TraversedDelCount;

                            if(ncVar.length() <= SHORT_DB_LENGTH)
                            {
                                ++TraversedShortDelCount;
                                ++impliedTICount;
                            }
                            else
                            {
                                impliedTICount += 2;
                            }

                            ++j; // past this DEL
                        }
                    }

                    ImpliedTICount += impliedTICount;
                }

                if(breakend.arm() == P_ARM && nextBreakend.arm() == Q_ARM)
                {
                    // if neither breakend points across the centromere then subtract it from the total range
                    // and also check for a reliable TI across this region as well
                    if(breakend.orientation() == 1 && nextBreakend.orientation() == -1 && !crossCentromereLinks.containsKey(chromosome))
                        TotalRange -= breakendsDistance;
                }
            }
        }
    }

    public static int findStartIndex(final List<SvBreakend> breakendList)
    {
        // first the first breakend in a cluster (on the chromosome), ignoring short TIs
        int startIndex = 0;
        for(; startIndex < breakendList.size() - 1; )
        {
            final SvBreakend lowerBreakend = breakendList.get(startIndex);
            final SvBreakend upperBreakend = breakendList.get(startIndex + 1);

            if(upperBreakend.position() - lowerBreakend.position() <= SHORT_TI_LENGTH
                    && lowerBreakend.orientation() == NEG_ORIENT && upperBreakend.orientation() == POS_ORIENT
                    && lowerBreakend.getLinkedPairs().stream().anyMatch(x -> upperBreakend.getLinkedPairs().contains(x)))
            {
                startIndex += 2;
                continue;
            }

            break;
        }

        return startIndex;
    }

    public static int findEndIndex(final List<SvBreakend> breakendList)
    {
        int endIndex = breakendList.size() - 1;

        for(; endIndex >= 1;)
        {
            final SvBreakend lowerBreakend = breakendList.get(endIndex - 1);
            final SvBreakend upperBreakend = breakendList.get(endIndex);

            if(upperBreakend.position() - lowerBreakend.position() <= SHORT_TI_LENGTH
                    && lowerBreakend.orientation() == NEG_ORIENT && upperBreakend.orientation() == POS_ORIENT
                    && lowerBreakend.getLinkedPairs().stream().anyMatch(x -> upperBreakend.getLinkedPairs().contains(x)))
            {
                endIndex -= 2;
                continue;
            }

            break;
        }

        return endIndex;
    }

}
