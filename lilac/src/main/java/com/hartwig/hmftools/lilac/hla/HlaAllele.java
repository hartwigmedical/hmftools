package com.hartwig.hmftools.lilac.hla;

import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_PREFIX;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class HlaAllele implements Comparable<HlaAllele>
{
    public final String Gene;
    public final String AlleleGroup;
    public final String Protein;
    public final String Synonymous;
    public final String SynonymousNonCoding;

    private final int mProteinNumber;
    private final HlaAllele mFourDigit;
    private final HlaAllele mGroup;
    private final int mHashCode;;
    private boolean mHasWildcards;

    public HlaAllele(
            final String gene, final String alleleGroup, final String protein, final String synonymous, final String synonymousNonCoding,
            final HlaAllele fourDigit, final HlaAllele group)
    {
        Gene = gene;
        AlleleGroup = alleleGroup;
        Protein = protein;
        Synonymous = synonymous;
        SynonymousNonCoding = synonymousNonCoding;

        mFourDigit = fourDigit;
        mGroup = group;

        mHashCode = toString().hashCode();
        mHasWildcards = false;

        int proteinNumber = 0;

        if(!protein.isEmpty())
        {
            try
            {
                proteinNumber = Integer.parseInt(protein);
            }
            catch(NumberFormatException e)
            {
                proteinNumber = Integer.parseInt(protein.substring(0, protein.length() - 1));
            }
        }

        mProteinNumber = proteinNumber;
    }

    public static HlaAllele fromString(final String line)
    {
        int starIndex = line.indexOf("*");
        String gene = line.substring(0, starIndex);
        String contigRemainder = line.substring(starIndex + 1);
        String[] contigSplit = contigRemainder.split(":");

        String alleleGroup = contigSplit[0];

        String protein = contigSplit.length < 2 ? "" : contigSplit[1];

        String synonymousCoding = contigSplit.length < 3 ? "" : contigSplit[2];
        String synonymousNonCoding = contigSplit.length < 4 ? "": contigSplit[3];

        HlaAllele group = !protein.isEmpty() ?
                new HlaAllele(gene, alleleGroup, "", "", "", null, null) : null;

        HlaAllele fourDigit = !synonymousCoding.isEmpty() ?
                new HlaAllele(gene, alleleGroup, protein, "", "", null, null) : null;

        return new HlaAllele(gene, alleleGroup, protein, synonymousCoding, synonymousNonCoding, fourDigit, group);
    }

    public String geneName() { return HLA_PREFIX + Gene; }

    public final HlaAllele asAlleleGroup()
    {
        return mGroup != null ? mGroup : this;
    }

    public final HlaAllele asFourDigit() { return mFourDigit != null ? mFourDigit : this; }

    public boolean hasWildcards() { return mHasWildcards; }
    public void setHasWildcard(boolean toggle) { mHasWildcards = toggle; }

    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if (!(other instanceof HlaAllele))
            return false;

        return hashCode() == other.hashCode();
    }

    public int hashCode() { return mHashCode; }
    public int proteinNumber() { return mProteinNumber; }

    @Override
    public int compareTo(final HlaAllele other)
    {
        if(hashCode() == other.hashCode())
            return 0;

        int geneCompare = Gene.compareTo(other.Gene);
        if(geneCompare != 0)
        {
            return geneCompare;
        }

        int groupCompare = AlleleGroup.compareTo(other.AlleleGroup);
        if(groupCompare != 0)
        {
            return groupCompare;
        }

        // force into integer terms rather than string (alphabetical) comparison
        if(proteinNumber() != other.proteinNumber())
        {
            return proteinNumber() > other.proteinNumber() ? 1 : -1;
        }

        int synonymousCodingCompare = Synonymous.compareTo(other.Synonymous);
        if(synonymousCodingCompare != 0)
        {
            return synonymousCodingCompare;
        }
        return SynonymousNonCoding.compareTo(other.SynonymousNonCoding);
    }

    public boolean matches(final HlaAllele other)
    {
        return other.toString().equals(toString());
    }

    public boolean matches(final String alleleStr)
    {
        return alleleStr.equals(toString());
    }

    public String toString()
    {
        CharSequence charSequence = Protein;
        if(charSequence.length() == 0)
        {
            return Gene + '*' + AlleleGroup;
        }
        charSequence = Synonymous;
        if(charSequence.length() == 0)
        {
            return Gene + '*' + AlleleGroup + ':' + Protein;
        }
        charSequence = SynonymousNonCoding;
        if(charSequence.length() == 0)
        {
            return Gene + '*' + AlleleGroup + ':' + Protein + ':' + Synonymous;
        }
        return Gene + '*' + AlleleGroup + ':' + Protein + ':' + Synonymous + ':' + SynonymousNonCoding;
    }

    public static String toString(final List<HlaAllele> alleles)
    {
        StringJoiner sj = new StringJoiner(", ");
        alleles.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public static List<HlaAllele> dedup(final List<HlaAllele> alleles)
    {
        // dedup but maintain ordering
        List<HlaAllele> newList = Lists.newArrayList();

        for(int i = 0; i < alleles.size(); ++i)
        {
            if(!newList.contains(alleles.get(i)))
                newList.add(alleles.get(i));
        }

        return newList;
    }

    public static Set<HlaAllele> findDuplicates(final List<HlaAllele> alleles)
    {
        Set<HlaAllele> duplicates = Sets.newHashSet();

        for(int i = 0; i < alleles.size() - 1; ++i)
        {
            for(int j = i + 1; j < alleles.size(); ++j)
            {
                if(alleles.get(i).matches(alleles.get(j)))
                {
                    duplicates.add(alleles.get(j));
                }
            }
        }

        return duplicates;
    }

}
