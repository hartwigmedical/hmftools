package com.hartwig.hmftools.lilac.nuc;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

public class NucleotideFragment
{
    private final String mId;
    private final Set<String> mGenes;
    private final List<Integer> mNucleotideLoci;
    private final List<Integer> mNucleotideQuality;
    private final List<String> mNucleotides;

    public final String getId() { return mId; }
    public final Set<String> getGenes() { return mGenes; }

    public final List<Integer> getNucleotideLoci() { return mNucleotideLoci; }
    public final List<Integer> getNucleotideQuality() { return mNucleotideQuality; }
    public final List<String> getNucleotides() { return mNucleotides; }

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

    public final boolean isEmpty() { return mNucleotideLoci.isEmpty(); }
    public final boolean isNotEmpty() { return !isEmpty(); }

    public int maxLoci() { return mNucleotideLoci.stream().mapToInt(x -> x).max().orElse(0); }

    public final boolean containsIndel()
    {
        return mNucleotides.stream().anyMatch(x -> x.equals(".") || x.length() > 1);
    }

    public final boolean containsNucleotide(int index)
    {
        return mNucleotideLoci.contains(index);
    }

    public final boolean containsAllNucleotides(final Set<Integer> indices)
    {
        return !indices.stream().anyMatch(x -> !containsNucleotide(x));
    }

    public final String nucleotides(final Set<Integer> indices)
    {
        StringJoiner sj = new StringJoiner("");
        indices.stream().forEach(x -> sj.add(nucleotide(x)));
        return sj.toString();
    }

    public final String nucleotide(int loci)
    {
        return mNucleotides.get(mNucleotideLoci.indexOf(loci));
    }

    public final NucleotideFragment qualityFilter(int minBaseQual)
    {
        // TODO: make static?
        final List<Integer> filteredNucleotideLoci = Lists.newArrayList();
        final List<Integer> filteredNucleotideQuality = Lists.newArrayList();
        final List<String> filteredNucleotides = Lists.newArrayList();

        // TODO
        return new NucleotideFragment(mId, mGenes, filteredNucleotideLoci, filteredNucleotideQuality, filteredNucleotides);
    }

    public AminoAcidFragment toAminoAcidFragment()
    {
        /*
                fun aminoAcid(index: Int): String {
            val first = nucleotide(index * 3)
            val second = nucleotide(index * 3 + 1)
            val third = nucleotide(index * 3 + 2)
            if (first == "." && second == "." && third == ".") {
                return "."
            }

            return Codons.aminoAcids(first + second + third)
        }

        val aminoAcidIndices = nucleotideLoci
                .filter { it % 3 == 0 }
                .filter { nucleotideLoci.contains(it + 1) && nucleotideLoci.contains(it + 2) }
                .map { it / 3 }

        val aminoAcids = aminoAcidIndices.map { aminoAcid(it) }.map { it }

        return AminoAcidFragment(id, genes, nucleotideLoci, nucleotideQuality, nucleotides, aminoAcidIndices, aminoAcids)

         */

        // TODO
        final List<Integer> aminoAcidLoci = Lists.newArrayList();
        final List<String> aminoAcids = Lists.newArrayList();

        return new AminoAcidFragment(mId, mGenes, mNucleotideLoci, mNucleotideQuality, mNucleotides, aminoAcidLoci, aminoAcids);
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
