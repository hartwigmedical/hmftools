package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacUtils.aboveMinQual;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.UNSET;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.read.Read;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

public class Fragment
{
    private final List<Read> mReads;
    private final String mReadGene; // mapped gene
    private final Set<String> mGenes; // other potentially applicable genes

    // initial nucleotide values
    private final List<Integer> mRawNucleotideLoci; // always in ascending order
    private final List<Integer> mRawNucleotideQuality;
    private final List<String> mRawNucleotides;

    // values which may be filtered
    private final List<Integer> mNucleotideLoci;
    private final List<Integer> mNucleotideQuality;
    private final List<String> mNucleotides;

    private boolean mIsQualFiltered;

    private int mAminoAcidConversionCount; // set true once nucleotides are converted into amino acids
    private final List<Integer> mAminoAcidLoci;
    private final List<String> mAminoAcids;

    private FragmentScope mScope;

    public Fragment(
            final Read read, final String readGene, final Set<String> genes, final List<Integer> nucleotideLoci,
            final List<Integer> qualities, final List<String> nucleotides)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mReads.add(read);

        mGenes = genes;
        mReadGene = readGene;

        mNucleotideLoci = Lists.newArrayList(nucleotideLoci);
        mNucleotideQuality = Lists.newArrayList(qualities);
        mNucleotides = Lists.newArrayList(nucleotides);

        mRawNucleotideLoci = Lists.newArrayList(nucleotideLoci);
        mRawNucleotideQuality = Lists.newArrayList(qualities);
        mRawNucleotides = Lists.newArrayList(nucleotides);

        mIsQualFiltered = false;

        mAminoAcidConversionCount = 0;
        mAminoAcidLoci = Lists.newArrayList();
        mAminoAcids = Lists.newArrayList();

        mScope = UNSET;

        if(mRawNucleotideLoci.size() != mRawNucleotides.size() || mRawNucleotideLoci.size() != mRawNucleotideQuality.size()
        || mNucleotides.size() != mNucleotideLoci.size() || mNucleotides.size() != mNucleotideQuality.size())
        {
            LL_LOGGER.warn("{} has differing raw loci counts(loci={} nuc={} quals={}) or nuc counts(loci={} nuc={} quals={})",
                    toString(), mNucleotideLoci.size(), mNucleotides.size(), mNucleotideQuality.size(),
                    mNucleotideLoci.size(), mNucleotides.size(), mNucleotideQuality.size());
        }
    }

    public void setAminoAcids(final List<Integer> aminoAcidLoci, final List<String> aminoAcids)
    {
        mAminoAcidLoci.addAll(aminoAcidLoci);
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

    public List<Integer> nucleotideLoci() { return mNucleotideLoci; }
    public List<Integer> nucleotideQuality() { return mNucleotideQuality; }
    public List<String> nucleotides() { return mNucleotides; }

    public List<Integer> rawNucleotideLoci() { return mRawNucleotideLoci; }
    public List<Integer> rawNucleotideQuality() { return mRawNucleotideQuality; }
    public List<String> rawNucleotides() { return mRawNucleotides; }

    public void addNucleotideInfo(int nucIndex, int locus, final String nucleotide, int quality)
    {
        if(nucIndex < 0 || nucIndex > mNucleotideLoci.size())
            return;

        // adds an extra base, quality and locus
        if(nucIndex == mNucleotideLoci.size() || mNucleotideLoci.get(nucIndex) != locus)
        {
            mNucleotideLoci.add(nucIndex, locus);
            mNucleotides.add(nucIndex, nucleotide);
            mNucleotideQuality.add(nucIndex, quality);

            if(!mIsQualFiltered)
            {
                // sets remain the same
                mRawNucleotideLoci.add(nucIndex, locus);
                mRawNucleotides.add(nucIndex, nucleotide);
                mRawNucleotideQuality.add(nucIndex, quality);
                return;
            }
        }

        // filtered indices have changed and are unlikely to match, so find the insertion point for this loci
        int rawIndex = findLociInsertionIndex(mRawNucleotideLoci, locus);

        if(rawIndex == mRawNucleotideLoci.size() || mRawNucleotideLoci.get(rawIndex) != locus)
        {
            mRawNucleotideLoci.add(rawIndex, locus);
            mRawNucleotides.add(rawIndex, nucleotide);
            mRawNucleotideQuality.add(rawIndex, quality);
        }
    }

    public void addNucleotideInfo(int locus, final String nucleotide, int quality)
    {
        int nucIndex = findLociInsertionIndex(mNucleotideLoci, locus);
        addNucleotideInfo(nucIndex, locus, nucleotide, quality);
    }

    public static int findLociInsertionIndex(final List<Integer> lociValues, final int locus)
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

    public boolean hasNucleotides() { return !mNucleotideLoci.isEmpty(); }

    public boolean containsIndel()
    {
        return mNucleotides.stream().anyMatch(x -> x.equals(".") || x.length() > 1);
    }

    public boolean containsNucleotide(int locus)
    {
        return mNucleotideLoci.contains(locus);
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
        int index = findLociIndex(mNucleotideLoci, locus);

        if(index < 0)
            return "";

        return mNucleotides.get(index);
    }

    public int minNucleotideLocus() { return !mNucleotideLoci.isEmpty() ? mNucleotideLoci.get(0) : -1; }
    public int maxNucleotideLocus() { return !mNucleotideLoci.isEmpty() ? mNucleotideLoci.get(mNucleotideLoci.size() - 1) : -1; }

    public List<Integer> aminoAcidLoci() { return mAminoAcidLoci; }
    public List<String> aminoAcids() { return mAminoAcids; }

    public int minAminoAcidLocus() { return !mAminoAcidLoci.isEmpty() ? mAminoAcidLoci.get(0) : -1; }
    public int maxAminoAcidLocus() { return !mAminoAcidLoci.isEmpty() ? mAminoAcidLoci.get(mAminoAcidLoci.size() - 1) : -1; }

    public void removeLowQualBases()
    {
        if(mIsQualFiltered)
            return;

        mIsQualFiltered = true;

        // cull any low-qual bases (they are retained in the raw bases)
        int index = 0;
        while(index < mNucleotideQuality.size())
        {
            if(aboveMinQual(mNucleotideQuality.get(index)))
            {
                ++index;
                continue;
            }

            mNucleotideQuality.remove(index);
            mNucleotideLoci.remove(index);
            mNucleotides.remove(index);
        }
    }

    public void buildAminoAcids()
    {
        ++mAminoAcidConversionCount;

        // build a amino-acid fragment from this nucleotide fragment, by creating amino acids for any complete codon
        mAminoAcids.clear();
        mAminoAcidLoci.clear();

        for(int i = 0; i < mNucleotideLoci.size(); ++i)
        {
            int locus = mNucleotideLoci.get(i);

            if(!isCodonMultiple(locus))
                continue;

            // since loci are ordered, can just check the next 2 expected do exist
            if(i >= mNucleotideLoci.size() - 2)
                break;

            if(mNucleotideLoci.get(i + 1) == locus + 1 && mNucleotideLoci.get(i + 2) == locus + 2)
            {
                int aaLocus = locus / 3;
                mAminoAcidLoci.add(aaLocus);
                mAminoAcids.add(formCodonAminoAcid(aaLocus));
            }
        }
    }

    public String formCodonAminoAcid(int locus)
    {
        return FragmentUtils.formCodonAminoAcid(locus, mNucleotideLoci, mNucleotides);
    }

    public boolean containsAminoAcids(final List<Integer> loci)
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

    public String aminoAcids(final List<Integer> loci)
    {
        StringJoiner sj = new StringJoiner("");
        loci.stream().forEach(x -> sj.add(aminoAcid(x)));
        return sj.toString();
    }

    public void filterOnLoci(final List<Integer> otherAminoAcidLoci)
    {
        int index = 0;
        while(index < mAminoAcidLoci.size())
        {
            int locus = mAminoAcidLoci.get(index);

            if(otherAminoAcidLoci.contains(locus))
            {
                ++index;
                continue;
            }

            mAminoAcidLoci.remove(index);
            mAminoAcids.remove(index);
        }
    }

    public String getLowQualNucleotide(int locus)
    {
        int index = mRawNucleotideLoci.indexOf(locus);
        return index >= 0 ? mRawNucleotides.get(index) : "";
    }

    public String getLowQualAminoAcid(int locus)
    {
        int startNucleotideLocus = locus * 3;

        int index = mRawNucleotideLoci.indexOf(startNucleotideLocus);

        if(index < 0 || index >= mRawNucleotideLoci.size() - 2)
            return "";

        if(mRawNucleotideLoci.get(index + 1) != startNucleotideLocus + 1 || mRawNucleotideLoci.get(index + 2) != startNucleotideLocus + 2)
            return "";

        return FragmentUtils.formCodonAminoAcid(locus, mRawNucleotideLoci, mRawNucleotides);
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
}
