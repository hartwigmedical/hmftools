package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacUtils.formRange;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import org.apache.commons.compress.utils.Lists;

public final class AminoAcidFragment extends NucleotideFragment
{
    private final List<Integer> mAminoAcidLoci;
    private final List<String> mAminoAcids;

    public AminoAcidFragment(
            final String id, final Set<String> genes, final List<Integer> nucleotideLoci, final List<Integer> nucleotideQuality, 
            final List<String> nucleotides, final List<Integer> aminoAcidLoci, final List<String> aminoAcids)
    {
        super(id, genes, nucleotideLoci, nucleotideQuality, nucleotides);
        mAminoAcidLoci = aminoAcidLoci;
        mAminoAcids = aminoAcids;
    }

    public List<Integer> getAminoAcidLoci() { return mAminoAcidLoci; }
    public List<String> getAminoAcids() { return mAminoAcids; }

    // cast to super class
    public static List<NucleotideFragment> nucFragments(final List<AminoAcidFragment> fragments)
    {
        return fragments.stream().collect(Collectors.toList());
    }

    public boolean containsAll(final List<Integer> loci)
    {
        return loci.stream().allMatch(x -> containsAminoAcid(x));
    }

    public boolean containsAminoAcid(int locus)
    {
        return mAminoAcidLoci.contains(locus);
    }

    public String aminoAcid(int locus)
    {
        int index = mAminoAcidLoci.indexOf(locus);

        if(index < 0)
            return "";

        return mAminoAcids.get(index);
    }

    public int maxAminoAcidLoci()
    {
        return mAminoAcidLoci.stream().mapToInt(x -> x).max().orElse(0);
    }

    // redundant
    public String aminoAcids(final List<Integer> loci)
    {
        StringJoiner sj = new StringJoiner("");
        loci.stream().forEach(x -> sj.add(aminoAcid(x)));
        return sj.toString();
    }

    public AminoAcidFragment intersectAminoAcidLoci(final List<Integer> otherAminoAcidLoci)
    {
        final List<Integer> intersectAminoAcidLoci = Lists.newArrayList();
        final List<String> intersectAminoAcids = Lists.newArrayList();

        boolean allPresent = true;

        for(int i = 0; i < mAminoAcidLoci.size(); ++i)
        {
            int locus = mAminoAcidLoci.get(i);

            if(otherAminoAcidLoci.contains(locus))
            {
                intersectAminoAcidLoci.add(locus);
                intersectAminoAcids.add(mAminoAcids.get(i));
            }
            else
            {
                allPresent = false;
            }
        }

        if(allPresent)
            return this;

        return new AminoAcidFragment(
                getId(), getGenes(), getNucleotideLoci(), getNucleotideQuality(), getNucleotides(),
                intersectAminoAcidLoci, intersectAminoAcids);
    }

    public String toString()
    {
        return String.format("%s genes(%s) lociRange(%d -> %d)",
                getId(), getGenes(), !mAminoAcidLoci.isEmpty() ? mAminoAcidLoci.get(0) : -1,
                !mAminoAcidLoci.isEmpty() ? mAminoAcidLoci.get(mAminoAcidLoci.size() - 1) : -1);
    }

    public boolean validate()
    {
        if(!super.validate())
            return false;

        // TODO - determine if valid or not, may still want to look back into nucleotide data and infer an AA
        //if(mAminoAcidLoci.isEmpty())
        //    return false;

        if(mAminoAcidLoci.size() != mAminoAcids.size())
            return false;

        for(int i = 0; i < mAminoAcidLoci.size() - 1; ++i)
        {
            if(mAminoAcidLoci.get(i) >= mAminoAcidLoci.get(i + 1))
                return false;
        }

        return true;
    }

    public static List<Integer> calcAminoAcidIndices(int nucStartIndex, int nucEndIndex)
    {
        int start = nucStartIndex / 3 + (nucStartIndex % 3 == 0 ? 0 : 1);
        int end = (nucEndIndex + 1) / 3 - 1;

        if(end <= start)
        {
            List<Integer> range = Lists.newArrayList();
            range.add(start);
            range.add(end);
            return range;
        }

        return formRange(start, end);
    }

}
