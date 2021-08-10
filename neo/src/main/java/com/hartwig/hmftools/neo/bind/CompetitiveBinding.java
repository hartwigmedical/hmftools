package com.hartwig.hmftools.neo.bind;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class CompetitiveBinding
{
    private final Map<String, Map<Integer, List<ScoreDistributionData>>> mAlleleRankDistributions;


    public CompetitiveBinding()
    {
        mAlleleRankDistributions = Maps.newHashMap();

        /*
        mScoreRankPercentileBuckets = Lists.newArrayListWithExpectedSize(16);

        double rankPercBucket = 0.00005;
        while(rankPercBucket < 1)
        {
            mScoreRankPercentileBuckets.add(rankPercBucket);
            rankPercBucket *= 2;
        }

         */
    }



}
