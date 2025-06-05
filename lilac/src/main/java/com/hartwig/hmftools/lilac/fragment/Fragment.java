package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacUtils.aboveMinQual;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.UNSET;

import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.read.Read;
import com.hartwig.hmftools.lilac.utils.AminoAcid;
import com.hartwig.hmftools.lilac.utils.Nucleotide;

import org.jetbrains.annotations.NotNull;

public class Fragment
{




    private final List<Read> mReads;
    private final String mReadGene; // mapped gene
    private final Set<String> mGenes; // other potentially applicable genes

    // initial nucleotide values, always in ascending order
    private final List<Nucleotide> mRawNucleotides;

    // values which may be filtered
    private final List<Nucleotide> mNucleotides;

    private boolean mIsQualFiltered;

    private int mAminoAcidConversionCount; // set true once nucleotides are converted into amino acids
    private final List<AminoAcid> mAminoAcids;

    private FragmentScope mScope;

    public Fragment(@NotNull final Read read, @NotNull final String readGene, @NotNull final Set<String> genes, @NotNull final List<Integer> nucleotideLoci, @NotNull final List<Byte> nucleotidesQualities, @NotNull final List<String> nucleotidesBases)
    {
        this(read, readGene, genes,
                IntStream.range(0, nucleotideLoci.size())
                        .mapToObj(i -> new Nucleotide(nucleotideLoci.get(i), nucleotidesQualities.get(i), nucleotidesBases.get(i)))
                        .collect(Collectors.toCollection(Lists::newArrayList)));
    }

    public Fragment(
            @NotNull final Read read, @NotNull final String readGene, @NotNull final Set<String> genes, @NotNull final List<Nucleotide> nucleotides)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mReads.add(read);

        mGenes = genes;
        mReadGene = readGene;

        mNucleotides = Lists.newArrayList(nucleotides);
        mRawNucleotides = Lists.newArrayList(nucleotides);

        mIsQualFiltered = false;

        mAminoAcids = Lists.newArrayList();

        mScope = UNSET;
    }

    // TODO: HERE
    public void setAminoAcids(final List<AminoAcid> aminoAcids)
    {
        mAminoAcids.addAll(aminoAcids);
    }

    public String id() { return mReads.get(0).Id; }

    public List<Read> reads() { return mReads; }

    public void addReads(final Fragment other)
    {
        other.reads().forEach(x -> addRead(x));
    }

    public void addRead(final Read read)
    {
        if(mReads.stream().anyMatch(x -> x.bamRecord().getFlags() == read.bamRecord().getFlags()))
            return;

        mReads.add(read);
    }

    public String readInfo() { return mReads.get(0).readInfo(); }
    public Set<String> genes() { return mGenes; }
    public void addGene(final String gene) { mGenes.add(gene); }
    public String readGene() { return mReadGene; }
    public boolean containsGene(final String gene) { return mGenes.stream().anyMatch(x -> x.equals(gene)); }

    public List<Nucleotide> nucleotides() { return mNucleotides; }
    public List<Nucleotide> rawNucleotides() { return mRawNucleotides; }

    public List<Integer> nucleotideLoci()
    {
        return Nucleotide.loci(mNucleotides);
    }

    public List<Integer> rawNucleotideLoci()
    {
        return Nucleotide.loci(mRawNucleotides);
    }

    public List<Integer> aminoAcidLoci()
    {
        return mAminoAcids.stream().map(AminoAcid::locus).toList();
    }

    public void addNucleotideInfo(int nucIndex, @NotNull final Nucleotide nucleotide)
    {
        if(nucIndex < 0 || nucIndex > mNucleotides.size())
            return;

        // adds an extra base, quality and locus
        if(nucIndex == mNucleotides.size() || mNucleotides.get(nucIndex).locus() != nucleotide.locus())
        {
            // TODO: random insertion into list.
            mNucleotides.add(nucIndex, nucleotide);

            if(!mIsQualFiltered)
            {
                // sets remain the same
                mRawNucleotides.add(nucIndex, nucleotide);
                return;
            }
        }

        // filtered indices have changed and are unlikely to match, so find the insertion point for this loci
        int rawIndex = findLociInsertionIndex(Nucleotide.loci(mRawNucleotides), nucleotide.locus());

        if(rawIndex == mRawNucleotides.size() || mRawNucleotides.get(rawIndex).locus() != nucleotide.locus())
            mRawNucleotides.add(rawIndex, nucleotide);
    }

    public void addNucleotideInfo(@NotNull final Nucleotide nucleotide)
    {
        int nucIndex = findLociInsertionIndex(Nucleotide.loci(mNucleotides), nucleotide.locus());
        addNucleotideInfo(nucIndex, nucleotide);
    }


    public void addNucleotideInfo(int nucIndex, int locus, @NotNull final String bases, byte quality)
    {
        Nucleotide nucleotide = new Nucleotide(locus, quality, bases);
        addNucleotideInfo(nucIndex, nucleotide);
    }

    public void addNucleotideInfo(int locus, @NotNull final String bases, byte quality)
    {
        Nucleotide nucleotide = new Nucleotide(locus, quality, bases);
        addNucleotideInfo(nucleotide);
    }

    public static int findLociInsertionIndex(@NotNull final List<Integer> lociValues, final int locus)
    {
        int index = Collections.binarySearch(lociValues, locus);

        if(index >= 0 && index < lociValues.size()) // check if matches
            return index;

        return -(index + 1); // as per insertion point conventions
    }

    public static int findLociIndex(final List<Integer> lociValues, final int locus)
    {
        int index = findLociInsertionIndex(lociValues, locus);

        if(index < 0 || index >= lociValues.size() || lociValues.get(index) != locus)
            return -1;

        return index;
    }

    public boolean hasNucleotides() { return !mNucleotides.isEmpty(); }

    public boolean containsIndel()
    {
        return mNucleotides.stream().map(Nucleotide::bases).anyMatch(x -> x.equals(".") || x.length() > 1);
    }

    public boolean containsNucleotide(int locus)
    {
        return mNucleotides.stream().mapToInt(Nucleotide::locus).anyMatch(x -> x == locus);
    }

    public boolean containsAllNucleotides(final List<Integer> loci)
    {
        return !loci.stream().anyMatch(x -> !containsNucleotide(x));
    }

    public String nucleotides(final List<Integer> indices)
    {
        StringJoiner sj = new StringJoiner("");
        indices.stream().forEach(x -> sj.add(nucleotide(x)));
        return sj.toString();
    }

    public String nucleotide(int locus)
    {
        int index = findLociIndex(Nucleotide.loci(mNucleotides), locus);

        if(index < 0)
            return "";

        return mNucleotides.get(index).bases();
    }

    public OptionalInt minNucleotideLocus() { return !mNucleotides.isEmpty() ? OptionalInt.of(mNucleotides.get(0).locus()) : OptionalInt.empty(); }
    public OptionalInt maxNucleotideLocus() { return !mNucleotides.isEmpty() ? OptionalInt.of(mNucleotides.get(mNucleotides.size() - 1).locus()) : OptionalInt.empty(); }

    public List<AminoAcid> aminoAcids() { return mAminoAcids; }

    public OptionalInt minAminoAcidLocus() { return !mAminoAcids.isEmpty() ? OptionalInt.of(mAminoAcids.get(0).locus()) : OptionalInt.empty(); }
    public OptionalInt maxAminoAcidLocus() { return !mAminoAcids.isEmpty() ? OptionalInt.of(mAminoAcids.get(mAminoAcids.size() - 1).locus()) : OptionalInt.empty(); }

    public void removeLowQualBases()
    {
        if(mIsQualFiltered)
            return;

        mIsQualFiltered = true;

        // cull any low-qual bases (they are retained in the raw bases)
        int index = 0;
        while(index < mNucleotides.size())
        {
            if(aboveMinQual(mNucleotides.get(index).qual()))
            {
                ++index;
                continue;
            }

            // TODO: removing from random position in a list.
            mNucleotides.remove(index);
        }
    }

    public void buildAminoAcids()
    {
        ++mAminoAcidConversionCount;

        // build an amino-acid fragment from this nucleotide fragment, by creating amino acids for any complete codon
        mAminoAcids.clear();

        for(int i = 0; i < mNucleotides.size(); ++i)
        {
            int locus = mNucleotides.get(i).locus();

            if(!isCodonMultiple(locus))
                continue;

            // since loci are ordered, can just check the next 2 expected do exist
            if(i >= mNucleotides.size() - 2)
                break;

            if(mNucleotides.get(i + 1).locus() == locus + 1 && mNucleotides.get(i + 2).locus() == locus + 2)
            {
                int aaLocus = locus / 3;
                mAminoAcids.add(new AminoAcid(aaLocus, formCodonAminoAcid(aaLocus)));
            }
        }
    }

    public String formCodonAminoAcid(int locus)
    {
        return FragmentUtils.formCodonAminoAcid(locus, mNucleotides);
    }

    public boolean containsAminoAcids(@NotNull final List<Integer> loci)
    {
        return loci.stream().allMatch(x -> containsAminoAcid(x));
    }

    public boolean containsAminoAcid(int locus)
    {
        return mAminoAcids.stream().mapToInt(AminoAcid::locus).anyMatch(x -> x == locus);
    }

    public String aminoAcid(int locus)
    {
        for(int i = 0; i < mAminoAcids.size(); i++)
        {
            var aminoAcid = mAminoAcids.get(i);
            if(aminoAcid.locus() == locus)
                return aminoAcid.acid();
        }

        return "";
    }

    public String aminoAcids(final List<Integer> loci)
    {
        StringJoiner sj = new StringJoiner("");
        loci.stream().forEach(x -> sj.add(aminoAcid(x)));
        return sj.toString();
    }

    public void filterOnLoci(@NotNull final List<Integer> otherAminoAcidLoci)
    {
        int index = 0;
        while(index < mAminoAcids.size())
        {
            int locus = mAminoAcids.get(index).locus();

            // TODO: has to search the whole list, multiple times
            if(otherAminoAcidLoci.contains(locus))
            {
                ++index;
                continue;
            }

            // TODO: removing from random position in a list
            mAminoAcids.remove(index);
        }
    }

    public String getLowQualNucleotide(int locus)
    {
        // TODO: skip searching the full list?
        for(int i = 0; i < mRawNucleotides.size(); i++)
        {
            Nucleotide nucleotide = mRawNucleotides.get(i);
            if(nucleotide.locus() == locus)
                return nucleotide.bases();
        }

        return "";
    }

    public String getLowQualAminoAcid(int locus)
    {
        int startNucleotideLocus = locus * 3;

        // TODO: fix having to stream.
        int index = mRawNucleotides.stream().map(Nucleotide::locus).toList().indexOf(startNucleotideLocus);

        if(index < 0 || index >= mRawNucleotides.size() - 2)
            return "";

        if(mRawNucleotides.get(index + 1).locus() != startNucleotideLocus + 1 || mRawNucleotides.get(index + 2).locus() != startNucleotideLocus + 2)
            return "";

        return FragmentUtils.formCodonAminoAcid(locus, mRawNucleotides);
    }

    public FragmentScope scope() { return mScope; }
    public void clearScope() { mScope = UNSET; }
    public boolean isScopeSet() { return mScope != UNSET; }

    public void setScope(FragmentScope scope) { setScope(scope, false); }

    public void setScope(FragmentScope scope, boolean override)
    {
        if(mScope != UNSET)
        {
            if(override)
            {
                if(scope != HLA_Y)
                {
                    LL_LOGGER.debug("frag({}: {}) overriding existing scope: {} -> {}", id(), readInfo(), mScope, scope);
                }

                mScope = scope;
            }
        }
        else
        {
            mScope = scope;
        }
    }

    @VisibleForTesting
    public int aminoAcidConversionCount() { return mAminoAcidConversionCount; }

    public String toString()
    {
        if(mAminoAcidConversionCount == 0)
        {
            return String.format("%s genes(%s) read(%s) nucRange(%d -> %d) qf(%s)",
                    id(), mGenes, readInfo(), minNucleotideLocus(), maxNucleotideLocus(), mIsQualFiltered);
        }

        return String.format("%s genes(%s) read(%s) nucRange(%d -> %d) qf(%s) aaRange(%d -> %d)",
                id(), mGenes, readInfo(), minNucleotideLocus(), maxNucleotideLocus(), mIsQualFiltered,
                minAminoAcidLocus(), maxAminoAcidLocus());
    }

    public boolean validate() { return FragmentUtils.validateFragment(this); }

    @VisibleForTesting
    public static int findLociIndexManual(final List<Integer> lociValues, final int locus)
    {
        int index = 0;

        while(index < lociValues.size())
        {
            if(locus < lociValues.get(index) || locus == lociValues.get(index))
                return index;

            if(lociValues.get(index) > locus)
                break;

            ++index;
        }

        return index;
    }

    // TODO: static analysis.
}
