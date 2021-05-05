package com.hartwig.hmftools.lilac.hla;

public class HlaAllele implements Comparable<HlaAllele>
{
    public final String Gene;
    public final String AlleleGroup;
    public final String Protein;
    public final String Synonymous;
    public final String SynonymousNonCoding;

    public HlaAllele(
            final String gene, final String alleleGroup, final String protein, final String synonymous, final String synonymousNonCoding)
    {
        Gene = gene;
        AlleleGroup = alleleGroup;
        Protein = protein;
        Synonymous = synonymous;
        SynonymousNonCoding = synonymousNonCoding;
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

    public final HlaAllele asSixDigit()
    {
        return new HlaAllele(Gene, AlleleGroup, Protein, Synonymous, "");
    }

    public final HlaAllele asFourDigit()
    {
        return new HlaAllele(Gene, AlleleGroup, Protein, "", "");
    }

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


    public boolean matches(final String alleleStr)
    {
        return alleleStr.equals(toString());
    }

}
