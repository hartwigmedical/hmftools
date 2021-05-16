package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Codons;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

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
        final List<Integer> filteredNucleotideLoci = Lists.newArrayList();
        final List<Integer> filteredNucleotideQuality = Lists.newArrayList();
        final List<String> filteredNucleotides = Lists.newArrayList();

        for(int i = 0; i < mNucleotideQuality.size(); ++i)
        {
            if(mNucleotideQuality.get(i) >= minBaseQual)
            {
                filteredNucleotideLoci.add(mNucleotideLoci.get(i));
                filteredNucleotideQuality.add(mNucleotideQuality.get(i));
                filteredNucleotides.add(mNucleotides.get(i));
            }
        }

        return new NucleotideFragment(mId, mGenes, filteredNucleotideLoci, filteredNucleotideQuality, filteredNucleotides);
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

    public final NucleotideFragment enrich(int loci, final String nucleotide, int quality)
    {
        final List<Integer> nucleotideLoci = Lists.newArrayList();
        nucleotideLoci.addAll(mNucleotideLoci);
        nucleotideLoci.add(loci);

        final List<Integer> nucleotideQuality = Lists.newArrayList();
        nucleotideQuality.addAll(mNucleotideQuality);
        nucleotideQuality.add(quality);

        final List<String> nucleotides = Lists.newArrayList();
        nucleotides.addAll(mNucleotides);
        nucleotides.add(nucleotide);

        return new NucleotideFragment(mId, mGenes, nucleotideLoci, nucleotideQuality, nucleotides);
    }

    public String toString()
    {
        return String.format("%s genes(%s) lociRange(%d -> %d)",
                mId, mGenes, !mNucleotideLoci.isEmpty() ? mNucleotideLoci.get(0) : -1,
                !mNucleotideLoci.isEmpty() ? mNucleotideLoci.get(mNucleotideLoci.size() - 1) : -1);
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
        final Set<String> genes = Sets.newHashSet();
        genes.addAll(frag1.getGenes());
        genes.addAll(frag2.getGenes());

        final List<Integer> nucleotideLoci = Lists.newArrayList();
        nucleotideLoci.addAll(frag1.getNucleotideLoci());
        nucleotideLoci.addAll(frag2.getNucleotideLoci());

        final List<Integer> nucleotideQuality = Lists.newArrayList();
        nucleotideQuality.addAll(frag1.getNucleotideQuality());
        nucleotideQuality.addAll(frag2.getNucleotideQuality());

        final List<String> nucleotides = Lists.newArrayList();
        nucleotides.addAll(frag1.getNucleotides());
        nucleotides.addAll(frag2.getNucleotides());

        return new NucleotideFragment(frag1.getId(), genes, nucleotideLoci, nucleotideQuality, nucleotides);
    }
}
