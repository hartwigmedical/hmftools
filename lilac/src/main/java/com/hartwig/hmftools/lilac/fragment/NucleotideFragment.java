package com.hartwig.hmftools.lilac.fragment;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Codons;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

public class NucleotideFragment
{
    private final String mId; // BAM read Id
    private final Set<String> mGenes;
    private final List<Integer> mNucleotideLoci;
    private final List<Integer> mNucleotideQuality;
    private final List<String> mNucleotides;

    public NucleotideFragment(
            final String id, final Set<String> genes, final List<Integer> nucleotideLoci,
            final List<Integer> nucleotideQuality, final List<String> nucleotides)
    {
        mId = id;
        mGenes = genes;
        mNucleotideLoci = nucleotideLoci;
        mNucleotideQuality = nucleotideQuality;
        mNucleotides = nucleotides;
    }

    public final String getId() { return mId; }
    public Set<String> getGenes() { return mGenes; }
    public boolean containsGene(final String gene) { return mGenes.stream().anyMatch(x -> x.equals(gene)); }

    public final List<Integer> getNucleotideLoci() { return mNucleotideLoci; }
    public final List<Integer> getNucleotideQuality() { return mNucleotideQuality; }
    public final List<String> getNucleotides() { return mNucleotides; }

    public final boolean isEmpty()
    {
        return mNucleotideLoci.isEmpty();
    }

    public final boolean isNotEmpty()
    {
        return !isEmpty();
    }

    public int maxLoci()
    {
        return mNucleotideLoci.stream().mapToInt(x -> x).max().orElse(0);
    }

    public final boolean containsIndel()
    {
        return mNucleotides.stream().anyMatch(x -> x.equals(".") || x.length() > 1);
    }

    public final boolean containsNucleotide(int index)
    {
        return mNucleotideLoci.contains(index);
    }

    public final boolean containsAllNucleotides(final List<Integer> indices)
    {
        return !indices.stream().anyMatch(x -> !containsNucleotide(x));
    }

    public final String nucleotides(final List<Integer> indices)
    {
        StringJoiner sj = new StringJoiner("");
        indices.stream().forEach(x -> sj.add(nucleotide(x)));
        return sj.toString();
    }

    public final String nucleotide(int locus)
    {
        int index = mNucleotideLoci.indexOf(locus);

        if(index < 0)
            return "";

        return mNucleotides.get(index);
    }

    public final NucleotideFragment qualityFilter(int minBaseQual)
    {
        final List<Integer> filteredIndices = Lists.newArrayList();
        boolean allPresent = true;

        for(int i = 0; i < mNucleotideQuality.size(); ++i)
        {
            if(mNucleotideQuality.get(i) >= minBaseQual)
            {
                filteredIndices.add(i);
            }
            else
            {
                allPresent = false;
            }
        }

        if(allPresent)
            return this;

        int filteredCount = filteredIndices.size();
        final List<Integer> filteredLoci = Lists.newArrayListWithExpectedSize(filteredCount);
        final List<Integer> filteredQuality = Lists.newArrayListWithExpectedSize(filteredCount);
        final List<String> filteredNucleotides = Lists.newArrayListWithExpectedSize(filteredCount);

        for(Integer index : filteredIndices)
        {
            filteredLoci.add(mNucleotideLoci.get(index));
            filteredQuality.add(mNucleotideQuality.get(index));
            filteredNucleotides.add(mNucleotides.get(index));
        }

        return new NucleotideFragment(mId, mGenes, filteredLoci, filteredQuality, filteredNucleotides);
    }

    public AminoAcidFragment toAminoAcidFragment()
    {
        final List<Integer> aminoAcidLoci = Lists.newArrayList();
        final List<String> aminoAcids = Lists.newArrayList();

        for(int i = 0; i < mNucleotideLoci.size(); ++i)
        {
            int locus = mNucleotideLoci.get(i);

            if((locus % 3) != 0)
            {
                continue;
            }

            if(mNucleotideLoci.contains(locus + 1) && mNucleotideLoci.contains(locus + 2))
            {
                int aaLocus = locus / 3;
                aminoAcidLoci.add(aaLocus);
                aminoAcids.add(formCodonAminoAcid(aaLocus));
            }
        }

        return new AminoAcidFragment(mId, mGenes, mNucleotideLoci, mNucleotideQuality, mNucleotides, aminoAcidLoci, aminoAcids);
    }

    private String formCodonAminoAcid(int index)
    {
        String first = nucleotide(index * 3);
        String second = nucleotide(index * 3 + 1);
        String third = nucleotide(index * 3 + 2);

        if(first == "." && second == "." && third == ".")
        {
            return ".";
        }

        return Codons.aminoAcids(first + second + third);
    }

    public void enrich(int locus, final String nucleotide, int quality)
    {
        // adds an extra base, quality and locus
        for(int i = 0; i < mNucleotideLoci.size(); ++i)
        {
            if(locus > mNucleotideLoci.get(i))
                continue;

            if(locus == mNucleotideLoci.get(i))
                break;

            // add in order
            mNucleotideLoci.add(i, locus);
            mNucleotideQuality.add(i, quality);
            mNucleotides.add(i, nucleotide);
            break;
        }
    }

    public String toString()
    {
        return String.format("%s genes(%s) lociRange(%d -> %d)",
                mId, mGenes, !mNucleotideLoci.isEmpty() ? mNucleotideLoci.get(0) : -1,
                !mNucleotideLoci.isEmpty() ? mNucleotideLoci.get(mNucleotideLoci.size() - 1) : -1);
    }

    public boolean validate()
    {
        if(mGenes.isEmpty() || mGenes.size() > 3)
            return false;

        // check loci are ordered and consistent with qualities and bases
        if(mNucleotides.isEmpty())
            return false;

        if(mNucleotides.size() != mNucleotideQuality.size() || mNucleotides.size() != mNucleotideLoci.size())
            return false;

        for(int i = 0; i < mNucleotideLoci.size() - 1; ++i)
        {
            if(mNucleotideLoci.get(i) >= mNucleotideLoci.get(i + 1))
                return false;
        }

        return true;
    }

    public static List<NucleotideFragment> reduceById(final List<NucleotideFragment> fragments)
    {
        // merge paired reads but keep any single reads as well
        final List<NucleotideFragment> mergedFragments = Lists.newArrayList(); // by readId

        Map<String,List<NucleotideFragment>> readGroupFrags = Maps.newHashMap();

        for(NucleotideFragment fragment : fragments)
        {
            List<NucleotideFragment> idFrags = readGroupFrags.get(fragment.getId());

            if(idFrags == null)
            {
                idFrags = Lists.newArrayList();
                readGroupFrags.put(fragment.getId(), idFrags);
            }

            idFrags.add(fragment);
        }

        for(Map.Entry<String,List<NucleotideFragment>> readGroup : readGroupFrags.entrySet())
        {
            List<NucleotideFragment> idFrags = readGroup.getValue();

            if(idFrags.size() == 1)
            {
                mergedFragments.add(idFrags.get(0));
            }
            else if(idFrags.size() == 2)
            {
                mergedFragments.add(NucleotideFragment.merge(idFrags.get(0), idFrags.get(1)));
            }
            else
            {
                NucleotideFragment mergedFrag = NucleotideFragment.merge(idFrags.get(0), idFrags.get(1));

                for(int i = 2; i < idFrags.size(); ++i)
                {
                    mergedFrag = NucleotideFragment.merge(mergedFrag, idFrags.get(i));
                }

                mergedFragments.add(mergedFrag);
            }
        }

        return mergedFragments;
    }

    public static List<Integer> expandIndices(final List<Integer> aminoAcidLoci)
    {
        List<Integer> nucleotideLoci = Lists.newArrayList();
        aminoAcidLoci.stream().forEach(x -> nucleotideLoci.addAll(Lists.newArrayList(3 * x, 3 * x + 1, 3 * x + 2)));
        return nucleotideLoci;
    }

    public static NucleotideFragment merge(final NucleotideFragment frag1, final NucleotideFragment frag2)
    {
        // merge frag 2 into 1 and ensure no repetition of loci
        frag2.getGenes().forEach(x -> frag1.getGenes().add(x));

        for(int index2 = 0; index2 < frag2.getNucleotideLoci().size(); ++index2)
        {
            int locus2 = frag2.getNucleotideLoci().get(index2);

            int index1 = 0;
            boolean add = true;
            for(; index1 < frag1.getNucleotideLoci().size(); ++index1)
            {
                int locus1 = frag1.getNucleotideLoci().get(index1);

                if(locus1 < locus2)
                    continue;

                if(locus1 == locus2)
                    add = false;

                break;
            }

            if(add)
            {
                frag1.getNucleotideLoci().add(index1, locus2);
                frag1.getNucleotides().add(index1, frag2.getNucleotides().get(index2));
                frag1.getNucleotideQuality().add(index1, frag2.getNucleotideQuality().get(index2));
            }
        }

        return frag1;
    }
}
