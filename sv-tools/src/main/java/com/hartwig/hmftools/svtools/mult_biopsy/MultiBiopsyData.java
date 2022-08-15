package com.hartwig.hmftools.svtools.mult_biopsy;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.purple.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

public class MultiBiopsyData
{
    public final String SampleId;

    public final int SvId;
    public final String ChrStart;
    public final String ChrEnd;
    public final long PosStart;
    public final long PosEnd;
    public final byte OrientStart;
    public final byte OrientEnd;
    public final StructuralVariantType Type;

    public final int ClusterId;
    public final int ClusterCount;
    public final String ResolvedType;

    private final List<MultiBiopsyData> mSharedMatches;
    private final List<MultiBiopsyData> mPartialMatches;

    private final List<String> mClusterReasons;

    public static final String MATCH_TYPE_PRIVATE = "Private";
    public static final String MATCH_TYPE_SHARED = "Shared";
    public static final String MATCH_TYPE_PARTIAL = "Partial";
    public static final String MATCH_TYPE_MIXED = "Mixed";

    public MultiBiopsyData(final String[] inputItems)
    {
        // SampleId,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,
        // ClusterId,ClusterCount,ClusterReason,ResolvedType,Ploidy,CNStart,CNEnd

        SampleId = inputItems[0];
        SvId = Integer.parseInt(inputItems[1]);

        String svType = inputItems[2];
        Type = svType.equals(NONE.toString()) ? SGL : StructuralVariantType.fromAttribute(svType);
        ChrStart = inputItems[3];
        PosStart = Long.parseLong(inputItems[4]);
        OrientStart = Byte.parseByte(inputItems[5]);
        ChrEnd = inputItems[6];
        PosEnd = Long.parseLong(inputItems[7]);
        OrientEnd = Byte.parseByte(inputItems[8]);

        ClusterId = Integer.parseInt(inputItems[9]);
        ClusterCount = Integer.parseInt(inputItems[10]);
        ResolvedType = inputItems[12];

        String clusterReasonsData = inputItems[11];
        mClusterReasons = Arrays.stream(clusterReasonsData.split(";")).collect(Collectors.toList());

        mSharedMatches = Lists.newArrayList();
        mPartialMatches = Lists.newArrayList();
    }

    public final String getMatchType()
    {
        if(!mSharedMatches.isEmpty() && !mPartialMatches.isEmpty())
            return MATCH_TYPE_MIXED;
        else if(!mSharedMatches.isEmpty())
            return MATCH_TYPE_SHARED;
        else if(!mPartialMatches.isEmpty())
            return MATCH_TYPE_PARTIAL;
        else
            return MATCH_TYPE_PRIVATE;

    }

    public void addMatchType(final String type, MultiBiopsyData other)
    {
        if(type == MATCH_TYPE_PARTIAL)
            mPartialMatches.add(other);
        else if(type == MATCH_TYPE_SHARED)
            mSharedMatches.add(other);
    }

    public long positionDiff(final MultiBiopsyData other)
    {
        return abs(PosStart - other.PosStart) + abs(PosEnd - other.PosEnd);
    }

    private static final int MAX_POS_DIFF = 20;

    public static boolean areMatched(String chr1, long pos1, byte orient1, String chr2, long pos2, byte orient2)
    {
        if (!chr1.equals(chr2) || orient1 != orient2)
            return false;

        long posDiff = abs(pos1 - pos2);
        return posDiff <= MAX_POS_DIFF;
    }

    public final List<MultiBiopsyData> getSharedMatches() { return mSharedMatches; }
    public final List<MultiBiopsyData> getPartialMatches() { return mPartialMatches; }
    public final List<String> getClusterReasons() { return mClusterReasons; }

    public void cullNonExactMatches()
    {
        // cull any partials or non-exact if exact are found
        if(mSharedMatches.isEmpty())
            return;

        mPartialMatches.clear();

        if(mSharedMatches.size() == 1)
            return;

        // only keep the best shared match
        MultiBiopsyData bestMatch = mSharedMatches.get(0);
        long minDiff = positionDiff(bestMatch);

        for(int i = 1; i < mSharedMatches.size(); ++i)
        {
            MultiBiopsyData other = mSharedMatches.get(i);
            long diff = positionDiff(other);

            if(diff < minDiff)
            {
                bestMatch = other;
                minDiff = diff;
            }
        }

        int index = 0;
        while(index < mSharedMatches.size())
        {
            if(mSharedMatches.get(index) != bestMatch)
            {
                MultiBiopsyData other = mSharedMatches.get(index);
                mSharedMatches.remove(index);
                other.getSharedMatches().remove(this);
                continue;
            }

            ++index;
        }
    }

    public String getClusterReasonForId(int otherId)
    {
        for(String clusterReasonId : mClusterReasons)
        {
            String[] crPair = clusterReasonId.split("_");
            if(crPair.length != 2)
                continue;

            if(Integer.parseInt(crPair[1]) == otherId)
                return crPair[0];
        }

        return "";
    }
}
