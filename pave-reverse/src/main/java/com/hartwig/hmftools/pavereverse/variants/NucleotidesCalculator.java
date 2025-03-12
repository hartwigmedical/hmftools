package com.hartwig.hmftools.pavereverse.variants;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;

public class NucleotidesCalculator
{
    private final AminoAcidSequence mAminoAcids;
    private final String mPrefix;
    private final String mSuffix;

    public NucleotidesCalculator(AminoAcidSequence aminoAcids, String prefix, String suffix)
    {
        Preconditions.checkArgument(prefix.length() < 3);
        Preconditions.checkArgument(suffix.length() < 3);
        Preconditions.checkArgument((prefix.length() + suffix.length()) <= (3 * aminoAcids.sequence().length()));
        mAminoAcids = aminoAcids;
        mPrefix = prefix;
        mSuffix = suffix;
    }

    public String anyBaseSequence()
    {
        // We only need this for situations where there is neither prefix nor suffix,
        // and we can only guarantee a non-null result with these preconditions.
        Preconditions.checkArgument(mPrefix.isEmpty());
        Preconditions.checkArgument(mSuffix.isEmpty());

        StringBuilder builder = new StringBuilder();
        for(int i = 0; i < mAminoAcids.length(); i++)
        {
            Set<String> allCodons = mAminoAcids.get(i).matchingTruncatedCodons("", "");
            // We want this to be deterministic, even if it's arbitrary
            String first = new TreeSet<>(allCodons).first();
            builder.append(first);
        }
        return builder.toString();
    }

    public Set<String> allPossibleBaseSequences()
    {
        List<Set<String>> candidateCodons = candidateAlternativeTruncatedCodons();
        // If there are no options for any of the codons, return the empty set.
        for(Set<String> codon : candidateCodons)
        {
            if(codon.isEmpty())
            {
                return new HashSet<>();
            }
        }
        return combinations(candidateAlternativeTruncatedCodons());
    }

    @VisibleForTesting
    List<Set<String>> candidateAlternativeTruncatedCodons()
    {
        final ArrayList<Set<String>> result = new ArrayList<>();
        for(int i = 0; i < mAminoAcids.length(); i++)
        {
            String requiredPrefix = prefixFilter(i);
            String requiredSuffix = suffixFilter(i);
            result.add(mAminoAcids.get(i).matchingTruncatedCodons(requiredPrefix, requiredSuffix));
        }
        return result;
    }

    private static Set<String> combinations(List<Set<String>> stringSets)
    {
        Set<String> leftSet = stringSets.get(0);
        if(stringSets.size() == 1)
        {
            return leftSet;
        }
        Set<String> remainderCombinations = combinations(stringSets.subList(1, stringSets.size()));
        Set<String> result = new HashSet<>();
        leftSet.forEach(l -> remainderCombinations.forEach(r -> result.add(l + r)));
        return result;
    }

    private String prefixFilter(int i)
    {
        if(i == 0)
        {
            return mPrefix;
        }
        return "";
    }

    private String suffixFilter(int i)
    {
        if(i == mAminoAcids.length() - 1)
        {
            return mSuffix;
        }
        return "";
    }
}
