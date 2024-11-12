package com.hartwig.hmftools.fastqtools.biomodalcollapse;

import static java.lang.Math.min;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.Nullable;

public class HairpinInfo
{
    private static final int KMER_MATCH_LENGTH = 8;
    private static final int MIN_KMER_MATCHES = 3;
    private static final int MIN_SUFFIX_MATCH_LENGTH = 4;

    public final int StartIndex;
    public final int MatchCount;
    public final int SuffixMatchLength;

    public HairpinInfo(final int startIndex, int matchCount, int suffixMatchLength)
    {
        StartIndex = startIndex;
        MatchCount = matchCount;
        SuffixMatchLength = suffixMatchLength;
    }

    @Nullable
    public static HairpinInfo findHairpin(final List<BaseQualPair> seq, final String hairpinSequence)
    {
        Map<Integer, Integer> matchCounts = Maps.newHashMap();
        int start = 0;
        int end = start + KMER_MATCH_LENGTH - 1;
        while(end < hairpinSequence.length())
        {
            int matchIndex = -1;
            for(int i = 0; i < seq.size() - (KMER_MATCH_LENGTH - 1); i++)
            {
                boolean allMatches = true;
                for(int j = 0; j < KMER_MATCH_LENGTH; j++)
                {
                    byte hairpinBase = (byte) hairpinSequence.charAt(start + j);
                    byte readBase = seq.get(i + j).Base;
                    if(hairpinBase != readBase)
                    {
                        allMatches = false;
                        break;
                    }
                }

                if(allMatches)
                {
                    matchIndex = i;
                    break;
                }
            }

            if(matchIndex < 0)
            {
                start++;
                end++;
                continue;
            }

            matchIndex -= start;
            matchCounts.put(matchIndex, matchCounts.getOrDefault(matchIndex, 0) + 1);

            start++;
            end++;
        }

        if(!matchCounts.isEmpty())
        {
            int bestIndex = -1;
            int bestCount = -1;
            for(Map.Entry<Integer, Integer> indexCountPair : matchCounts.entrySet())
            {
                int index = indexCountPair.getKey();
                int count = indexCountPair.getValue();
                if(count > bestCount)
                {
                    bestIndex = index;
                    bestCount = count;
                }
            }

            if(bestCount >= MIN_KMER_MATCHES)
            {
                return new HairpinInfo(bestIndex, bestCount, -1);
            }
        }

        // look for exact matches at end of read
        int maxHairpinPrefixLength = min(hairpinSequence.length(), KMER_MATCH_LENGTH + MIN_KMER_MATCHES - 2);
        for(int hairpinPrefixLength = maxHairpinPrefixLength; hairpinPrefixLength >= MIN_SUFFIX_MATCH_LENGTH; hairpinPrefixLength--)
        {
            start = seq.size() - hairpinPrefixLength;
            if(start < 0)
            {
                continue;
            }

            boolean all_matches = true;
            for(int i = 0; i < hairpinPrefixLength; i++)
            {
                byte readBase = seq.get(start + i).Base;
                byte hairpinBase = (byte) hairpinSequence.charAt(i);
                if(readBase != hairpinBase)
                {
                    all_matches = false;
                    break;
                }
            }

            if(all_matches)
            {
                return new HairpinInfo(start, -1, hairpinPrefixLength);
            }
        }

        return null;
    }
}
