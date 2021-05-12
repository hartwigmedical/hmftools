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

    public HlaAllele(
            final String gene, final String alleleGroup, final String protein, final String synonymous, final String synonymousNonCoding)
    {
        Gene = gene;
        AlleleGroup = alleleGroup;
        Protein = protein;
        Synonymous = synonymous;
        SynonymousNonCoding = synonymousNonCoding;

        if(!synonymous.isEmpty())
            mFourDigit = new HlaAllele(gene, alleleGroup, protein, "", "");
        else
            mFourDigit = null;
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

        return new HlaAllele(gene, alleleGroup, protein, synonymousCoding, synonymousNonCoding);
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

    public static String toString(final List<HlaAllele> allees)
    {
        StringJoiner sj = new StringJoiner(", ");
        allees.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public final HlaAllele asSixDigit()
    {
        return new HlaAllele(Gene, AlleleGroup, Protein, Synonymous, "");
    }

    public final HlaAllele asFourDigit() { return mFourDigit != null ? mFourDigit : this; }

    public final HlaAllele asAlleleGroup()
    {
        return new HlaAllele(Gene, AlleleGroup, "", "", "");
    }

    @Override
    public int compareTo(final HlaAllele other)
    {
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

    public static boolean matches(final List<HlaAllele> list1, final List<HlaAllele> list2)
    {
        if(list1.size() != list2.size())
            return false;

        return list1.stream().allMatch(x -> contains(list2, x));
    }

    public static boolean contains(final List<HlaAllele> list, final HlaAllele allele)
    {
        return list.stream().anyMatch(x -> x.matches(allele));
    }

    public static List<HlaAllele> dedup(final List<HlaAllele> alleles)
    {
        List<HlaAllele> newList = Lists.newArrayList();

        for(int i = 0; i < alleles.size(); ++i)
        {
            if(!contains(newList, alleles.get(i)))
                newList.add(alleles.get(i));
        }

        return newList;
    }

}
