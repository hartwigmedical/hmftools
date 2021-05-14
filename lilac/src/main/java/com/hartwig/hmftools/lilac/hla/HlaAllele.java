package com.hartwig.hmftools.lilac.hla;

import static java.lang.Math.min;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public class HlaAllele implements Comparable<HlaAllele>
{
    public final String Gene;
    public final String AlleleGroup;
    public final String Protein;
    public final String Synonymous;
    public final String SynonymousNonCoding;

    private final HlaAllele mFourDigit;
    private final HlaAllele mGroup;
    private final int mHashCode;;

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

    public boolean isGroup() { return Protein.isEmpty(); }
    public boolean isFourDigit() { return Synonymous.isEmpty(); }

    public final HlaAllele asAlleleGroup()
    {
        return mGroup != null ? mGroup : this;
    }

    public final HlaAllele asFourDigit() { return mFourDigit != null ? mFourDigit : this; }

    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if (!(other instanceof HlaAllele))
            return false;

        return hashCode() == other.hashCode();
    }

    public int hashCode() { return mHashCode; }

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
        int proteinCompare = Protein.compareTo(other.Protein);
        if(proteinCompare != 0)
        {
            return proteinCompare;
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

    public static List<HlaAllele> takeN(final List<HlaAllele> list, int n)
    {
        List<HlaAllele> newList = Lists.newArrayList();

        for(int i = 0; i < min(list.size(), n); ++i)
        {
            newList.add(list.get(i));
        }

        return newList;
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

    public static boolean contains(final List<HlaAllele> list, final HlaAllele allele)
    {
        return list.contains(allele);
    }

    public static List<HlaAllele> dedup(final List<HlaAllele> alleles)
    {
        // dedep but maintain ordering
        List<HlaAllele> newList = Lists.newArrayList();

        for(int i = 0; i < alleles.size(); ++i)
        {
            if(!newList.contains(alleles.get(i)))
                newList.add(alleles.get(i));
        }

        return newList;
    }

}
