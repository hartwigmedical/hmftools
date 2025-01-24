package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.common.codon.AminoAcids.AMINO_ACID_TO_CODON_MAP;

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
    private final char[] mAminoAcids;

    @NotNull
    private final String prefix;

    @NotNull
    private final String suffix;

    public NucleotidesCalculator(@NotNull final String aminoAcids, @NotNull final String prefix, @NotNull final String suffix)
    {
        Preconditions.checkArgument(Checks.isValidProtein(aminoAcids));
        Preconditions.checkArgument(prefix.length() < 3);
        Preconditions.checkArgument(suffix.length() < 3);
        Preconditions.checkArgument((prefix.length() + suffix.length()) <= (3 * aminoAcids.length()));
        mAminoAcids = aminoAcids.toCharArray();
        this.prefix = prefix;
        this.suffix = suffix;
    }

    @NotNull
    public Set<String> possibilities()
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
        for(int i = 0; i < mAminoAcids.length; i++)
        {
            result.add(matchingTruncatedCodons(i));
        }
        return result;
    }

    private Set<String> matchingTruncatedCodons(int index)
    {
        String requiredPrefix = prefixFilter(index);
        String requiredSuffix = suffixFilter(index);
        Set<String> result = new HashSet<>();
        String aminoAcid = mAminoAcids[index] + "";
        AMINO_ACID_TO_CODON_MAP.get(aminoAcid).forEach(codon ->
        {
            if(codon.startsWith(requiredPrefix) && codon.endsWith(requiredSuffix))
            {
                String prefixTruncatedCodon = codon.substring(requiredPrefix.length());
                String truncatedCodon = prefixTruncatedCodon.substring(0, prefixTruncatedCodon.length() - requiredSuffix.length());
                if(!truncatedCodon.isEmpty())
                {
                    result.add(truncatedCodon);
                }
            }
        });
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
        if(i == mAminoAcids.length - 1)
        {
            return suffix;
        }
        return "";
    }
}
