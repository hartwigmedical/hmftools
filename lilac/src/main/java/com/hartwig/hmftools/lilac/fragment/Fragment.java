package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.UNSET;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.validateLociBases;

import com.google.common.collect.Lists;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class Fragment
{
    private final String mId; // BAM read Id
    private final String mReadInfo; // BAM read chr, start, end, gene, cigar

    private final String mReadGene; // mapped gene
    private final Set<String> mGenes;
    private final List<Integer> mNucleotideLoci;
    private final List<Integer> mNucleotideQuality;
    private final List<String> mNucleotides;

    private boolean mIsQualFiltered;
    private final List<Integer> mRawNucleotideLoci;
    private final List<Integer> mRawNucleotideQuality;
    private final List<String> mRawNucleotides;

    private int mAminoAcidConversionCount;
    private final List<Integer> mAminoAcidLoci;
    private final List<String> mAminoAcids;

    private FragmentScope mScope;

    public Fragment(
            final String id, final String readInfo, final String readGene, final Set<String> genes, final List<Integer> nucleotideLoci,
            final List<Integer> qualities, final List<String> nucleotides)
    {
        mId = id;
        mReadInfo = readInfo;

        mGenes = genes;
        mReadGene = readGene;

        mNucleotideLoci = nucleotideLoci.stream().collect(Collectors.toList());
        mNucleotideQuality = qualities.stream().collect(Collectors.toList());
        mNucleotides = nucleotides.stream().collect(Collectors.toList());

        // independent copies
        mRawNucleotideLoci = nucleotideLoci.stream().collect(Collectors.toList());
        mRawNucleotideQuality = qualities.stream().collect(Collectors.toList());
        mRawNucleotides = nucleotides.stream().collect(Collectors.toList());

        mIsQualFiltered = false;

        mAminoAcidConversionCount = 0;
        mAminoAcidLoci = Lists.newArrayList();
        mAminoAcids = Lists.newArrayList();

        mScope = UNSET;
    }

    public void setAminoAcids(final List<Integer> aminoAcidLoci, final List<String> aminoAcids)
    {
        mAminoAcidLoci.addAll(aminoAcidLoci);
        mAminoAcids.addAll(aminoAcids);
    }

    public String id() { return mId; }
    public String readInfo() { return mReadInfo; }
    public Set<String> getGenes() { return mGenes; }
    public String readGene() { return mReadGene; }
    public boolean containsGene(final String gene) { return mGenes.stream().anyMatch(x -> x.equals(gene)); }

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
                    LL_LOGGER.debug("frag({}: {}) overriding existing scope: {} -> {}", mId, mReadInfo, mScope, scope);
                }

                mScope = scope;
            }
        }
        else
        {
            mScope = scope;
        }
    }

    public List<Integer> getNucleotideLoci() { return mNucleotideLoci; }
    public List<Integer> getNucleotideQuality() { return mNucleotideQuality; }
    public List<String> getNucleotides() { return mNucleotides; }

    public List<Integer> getRawNucleotideLoci() { return mRawNucleotideLoci; }
    public List<Integer> getRawNucleotideQuality() { return mRawNucleotideQuality; }
    public List<String> getRawNucleotides() { return mRawNucleotides; }

    public void addNucleotide(int index, int locus, String bases, int quality)
    {
        mNucleotideLoci.add(index, locus);
        mNucleotides.add(index, bases);
        mNucleotideQuality.add(index, quality);

        if(mRawNucleotideLoci.size() == mNucleotideLoci.size())
        {
            mRawNucleotideLoci.add(index, locus);
            mRawNucleotides.add(index, bases);
            mRawNucleotideQuality.add(index, quality);
        }
        else
        {
            int rawIndex = 0;

            while(rawIndex < mRawNucleotideLoci.size())
            {
                if(mRawNucleotideLoci.get(rawIndex) > locus)
                    break;

                ++rawIndex;
            }

            mRawNucleotideLoci.add(rawIndex, locus);
            mRawNucleotides.add(rawIndex, bases);
            mRawNucleotideQuality.add(rawIndex, quality);
        }
    }

    public boolean hasNucleotides() { return !mNucleotideLoci.isEmpty(); }

    public boolean containsIndel()
    {
        return mNucleotides.stream().anyMatch(x -> x.equals(".") || x.length() > 1);
    }

    public boolean containsNucleotide(int index)
    {
        return mNucleotideLoci.contains(index);
    }

    public boolean containsAllNucleotides(final List<Integer> indices)
    {
        return !indices.stream().anyMatch(x -> !containsNucleotide(x));
    }

    public String nucleotides(final List<Integer> indices)
    {
        StringJoiner sj = new StringJoiner("");
        indices.stream().forEach(x -> sj.add(nucleotide(x)));
        return sj.toString();
    }

    public String nucleotide(int locus)
    {
        int index = mNucleotideLoci.indexOf(locus);

        if(index < 0)
            return "";

        return mNucleotides.get(index);
    }

    public int minNucleotideLocus() { return !mNucleotideLoci.isEmpty() ? mNucleotideLoci.get(0) : -1; }
    public int maxNucleotideLocus() { return !mNucleotideLoci.isEmpty() ? mNucleotideLoci.get(mNucleotideLoci.size() - 1) : -1; }
    public int maxLoci() { return maxNucleotideLocus(); }

    public List<Integer> getAminoAcidLoci() { return mAminoAcidLoci; }
    public List<String> getAminoAcids() { return mAminoAcids; }

    public int minAminoAcidLocus() { return !mAminoAcidLoci.isEmpty() ? mAminoAcidLoci.get(0) : -1; }
    public int maxAminoAcidLocus() { return !mAminoAcidLoci.isEmpty() ? mAminoAcidLoci.get(mAminoAcidLoci.size() - 1) : -1; }

    public void qualityFilter(int minBaseQual)
    {
        if(mIsQualFiltered)
            return;

        mIsQualFiltered = true;

        // cull any low-qual bases (they are retained in the raw bases)
        final List<Integer> filteredIndices = Lists.newArrayList();

        int index = 0;
        while(index < mNucleotideQuality.size())
        {
            if(mNucleotideQuality.get(index) >= minBaseQual)
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

    public void enrich(int locus, final String nucleotide, int quality)
    {
        // adds an extra base, quality and locus
        for(int i = 0; i < mNucleotideLoci.size(); ++i)
        {
            if(locus > mNucleotideLoci.get(i))
                continue;

            if(locus == mNucleotideLoci.get(i))
                break;

            addNucleotide(i, locus, nucleotide, quality);
            break;
        }
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

    public String toString()
    {
        if(mAminoAcidConversionCount == 0)
        {
            return String.format("%s genes(%s) read(%s) nucRange(%d -> %d) qf(%s)",
                    mId, mGenes, readInfo(), minNucleotideLocus(), maxNucleotideLocus(), mIsQualFiltered);
        }

        return String.format("%s genes(%s) read(%s) nucRange(%d -> %d) qf(%s) aaRange(%d -> %d)",
                mId, mGenes, readInfo(), minNucleotideLocus(), maxNucleotideLocus(), mIsQualFiltered,
                minAminoAcidLocus(), maxAminoAcidLocus());
    }

    public boolean validate()
    {
        if(mGenes.isEmpty() || mGenes.size() > 3)
        {
            LL_LOGGER.warn("{} {} has no genes", mId, mReadInfo);
            return false;
        }

        // check loci are ordered and consistent with qualities and bases
        if(mNucleotides.isEmpty() || mRawNucleotides.isEmpty())
        {
            LL_LOGGER.warn("{} {} has no bases", mId, mReadInfo);
            return false;
        }

        if(mNucleotides.size() != mNucleotideQuality.size())
        {
            LL_LOGGER.warn("{} {} inconsistent bases loci({}) quals({})",
                    mId, mReadInfo, mNucleotides.size(), mNucleotideQuality.size());
            return false;
        }

        if(!validateLociBases(mId, mNucleotideLoci, mNucleotides))
            return false;

        if(!validateLociBases(mId, mRawNucleotideLoci, mRawNucleotides))
            return false;

        if(mAminoAcidConversionCount > 0)
        {
            if(!validateLociBases(mId, mAminoAcidLoci, mAminoAcids))
                return false;

            for(int i = 0; i < mAminoAcidLoci.size(); ++i)
            {
                if(mAminoAcids.get(i).isEmpty())
                {
                    LL_LOGGER.warn("{} {} empty amino-acid, index({})", mId, mReadInfo, i);
                    return false;
                }
            }

            if(mAminoAcidConversionCount > 1)
            {
                LL_LOGGER.warn("{} {} amino-acid conversion repeated({})", mId, mReadInfo, mAminoAcidConversionCount);
                return false;
            }
        }

        return true;
    }

}
