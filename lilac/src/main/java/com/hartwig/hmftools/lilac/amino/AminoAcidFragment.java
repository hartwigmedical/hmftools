package com.hartwig.hmftools.lilac.amino;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.hartwig.hmftools.lilac.nuc.NucleotideFragment;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

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

    public final boolean containsAll(final List<Integer> indices)
    {
        return !indices.stream().noneMatch(x -> mAminoAcidLoci.contains(x));
    }

    public final boolean containsAminoAcid(int index)
    {
        return mAminoAcidLoci.contains(index);
    }

    public final String aminoAcid(int loci)
    {
        return mAminoAcids.get(mAminoAcidLoci.indexOf(loci));
    }

    @NotNull
    public final String aminoAcids(final List<Integer> indices)
    {
        StringJoiner sj = new StringJoiner("");
        indices.stream().forEach(x -> sj.add(aminoAcid(x)));
        return sj.toString();
    }

    public final List<Integer> aminoAcidLoci() { return mAminoAcidLoci; }

    public final AminoAcidFragment intersectAminoAcidLoci(final List<Integer> aminoAcidLoci)
    {
        /*
                val filteredIndexes = aminoAcidLoci
                .mapIndexed { index: Int, loci: Int -> Pair(index, loci) }
                .filter { (_, loci) -> loci in otherAminoAcidLoci }
                .map { (index, _) -> index }

        val filteredAminoAcidLoci = filteredIndexes.map { aminoAcidLoci[it] }
        val filteredAminoAcids = filteredIndexes.map { aminoAcids[it] }

        return AminoAcidFragment(id, genes, nucleotideLoci, nucleotideQuality, nucleotides, filteredAminoAcidLoci, filteredAminoAcids)

         */

        // TODO

        final List<Integer> intersectAminoAcidLoci = Lists.newArrayList();
        final List<String> intersectAminoAcids = Lists.newArrayList();


        return new AminoAcidFragment(
                getId(), getGenes(), getNucleotideLoci(), getNucleotideQuality(), getNucleotides(),
                intersectAminoAcidLoci, intersectAminoAcids);
    }
}
