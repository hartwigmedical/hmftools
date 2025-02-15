package com.hartwig.hmftools.pave.transval;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public class NucleotidesCalculator
{
    @NotNull
    private final AminoAcidSequence mAminoAcids;

    @NotNull
    private final String prefix;

    @NotNull
    private final String suffix;

    NucleotidesCalculator(@NotNull final AminoAcidSequence aminoAcids, @NotNull final String prefix, @NotNull final String suffix)
    {
        Preconditions.checkArgument(prefix.length() < 3);
        Preconditions.checkArgument(suffix.length() < 3);
        Preconditions.checkArgument((prefix.length() + suffix.length()) <= (3 * aminoAcids.sequence().length()));
        mAminoAcids = aminoAcids;
        this.prefix = prefix;
        this.suffix = suffix;
    }

    @NotNull
    public String anyBaseSequence()
    {
        // We only need this for situations where there is neither prefix nor suffix,
        // and we can only guarantee a non-null result with these preconditions.
        Preconditions.checkArgument(prefix.isEmpty());
        Preconditions.checkArgument(suffix.isEmpty());

        StringBuilder builder = new StringBuilder();
        for(int i = 0; i < mAminoAcids.length(); i++)
        {
            Set<String> allCodons = mAminoAcids.get(i).matchingTruncatedCodons("","");
            builder.append(allCodons.iterator().next());
        }
        return builder.toString();
    }

    @NotNull
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
    public List<Set<String>> candidateAlternativeTruncatedCodons()
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
            return prefix;
        }
        return "";
    }

    private String suffixFilter(int i)
    {
        if(i == mAminoAcids.length() - 1)
        {
            return suffix;
        }
        return "";
    }
}
