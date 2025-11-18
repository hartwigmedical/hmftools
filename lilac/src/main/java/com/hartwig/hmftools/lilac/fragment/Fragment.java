package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.UNSET;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.evidence.AminoAcid;
import com.hartwig.hmftools.lilac.evidence.Nucleotide;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.read.Read;

public class Fragment
{
    private final List<Read> mReads;
    private final HlaGene mReadGene; // mapped gene
    private final Set<HlaGene> mGenes; // other potentially applicable genes

    // initial nucleotide values
    private final NavigableMap<Integer, Nucleotide> mRawNucleotidesByLoci;

    // values which may be filtered
    private final NavigableMap<Integer, Nucleotide> mNucleotidesByLoci;

    private boolean mIsQualFiltered;

    private int mAminoAcidConversionCount; // set true once nucleotides are converted into amino acids
    private final NavigableMap<Integer, AminoAcid> mAminoAcidsByLoci;

    private FragmentScope mScope;

    public static Fragment createFromQuals(final Read read, final HlaGene readGene, final Set<HlaGene> genes,
            final List<Integer> nucleotideLoci, final List<Byte> nucleotidesQualities, final List<String> nucleotidesBases)
    {
        return new Fragment(read, readGene, genes,
                IntStream.range(0, nucleotideLoci.size())
                        .mapToObj(i -> Nucleotide.create(nucleotideLoci.get(i), nucleotidesQualities.get(i), nucleotidesBases.get(i)))
                        .toList());
    }

    public static Fragment createFromIsLowQuals(final Read read, final HlaGene readGene, final Set<HlaGene> genes,
            final List<Integer> nucleotideLoci, final List<Boolean> isLowQuals, final List<String> nucleotidesBases)
    {
        return new Fragment(read, readGene, genes,
                IntStream.range(0, nucleotideLoci.size())
                        .mapToObj(i -> Nucleotide.create(nucleotideLoci.get(i), isLowQuals.get(i), nucleotidesBases.get(i)))
                        .toList());
    }

    public Fragment(
            final Read read, final HlaGene readGene, final Set<HlaGene> genes, final Iterable<Nucleotide> nucleotides)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mReads.add(read);

        mGenes = genes;
        mReadGene = readGene;

        mNucleotidesByLoci = Maps.newTreeMap();
        mRawNucleotidesByLoci = Maps.newTreeMap();
        for(Nucleotide nucleotide : nucleotides)
        {
            mNucleotidesByLoci.put(nucleotide.locus(), nucleotide);
            mRawNucleotidesByLoci.put(nucleotide.locus(), nucleotide);
        }

        mIsQualFiltered = false;

        mAminoAcidsByLoci = Maps.newTreeMap();

        mScope = UNSET;
    }

    public String id() { return mReads.get(0).Id; }

    public List<Read> reads() { return mReads; }

    public void addReads(final Fragment other) { other.reads().forEach(this::addRead); }

    public void addRead(final Read read)
    {
        if(mReads.stream().anyMatch(x -> x.bamRecord().getFlags() == read.bamRecord().getFlags()))
            return;

        mReads.add(read);
    }

    public String readInfo() { return mReads.get(0).readInfo(); }
    public Set<HlaGene> genes() { return mGenes; }
    public void addGene(final HlaGene gene) { mGenes.add(gene); }
    public HlaGene readGene() { return mReadGene; }
    public boolean containsGene(final HlaGene gene) { return mGenes.contains(gene); }

    public NavigableMap<Integer, Nucleotide> nucleotidesByLoci() { return mNucleotidesByLoci; }
    public NavigableMap<Integer, Nucleotide> rawNucleotidesByLoci() { return mRawNucleotidesByLoci; }

    public void addNucleotide(final Nucleotide nucleotide)
    {
        // adds an extra base, quality and locus
        if(!mNucleotidesByLoci.containsKey(nucleotide.locus()))
        {
            mNucleotidesByLoci.put(nucleotide.locus(), nucleotide);
            if(!mIsQualFiltered)
            {
                // sets remain the same
                mRawNucleotidesByLoci.put(nucleotide.locus(), nucleotide);
                return;
            }
        }

        if(!mRawNucleotidesByLoci.containsKey(nucleotide.locus()))
            mRawNucleotidesByLoci.put(nucleotide.locus(), nucleotide);
    }

    public void addNucleotide(int locus, final String bases, byte quality)
    {
        Nucleotide nucleotide = Nucleotide.create(locus, quality, bases);
        addNucleotide(nucleotide);
    }

    public void addHighQualNucleotide(int locus, final String bases)
    {
        Nucleotide nucleotide = Nucleotide.createHighQual(locus, bases);
        addNucleotide(nucleotide);
    }

    public boolean hasNucleotides() { return !mNucleotidesByLoci.isEmpty(); }

    public boolean containsIndel()
    {
        return mNucleotidesByLoci.values().stream().map(Nucleotide::bases).anyMatch(x -> ".".equals(x) || x.length() > 1);
    }

    public boolean containsNucleotideLocus(int locus)
    {
        return mNucleotidesByLoci.containsKey(locus);
    }

    public boolean containsAllNucleotideLoci(final Collection<Integer> loci)
    {
        return loci.stream().allMatch(this::containsNucleotideLocus);
    }

    public String nucleotides(final Iterable<Integer> loci)
    {
        StringJoiner sj = new StringJoiner("");
        loci.forEach(x -> sj.add(nucleotide(x)));
        return sj.toString();
    }

    public String nucleotide(int locus)
    {
        Nucleotide nuc = mNucleotidesByLoci.get(locus);
        if(nuc == null)
            return "";

        return nuc.bases();
    }

    public int minNucleotideLocus() { return mNucleotidesByLoci.isEmpty() ? -1 : mNucleotidesByLoci.firstKey(); }
    public int maxNucleotideLocus() { return mNucleotidesByLoci.isEmpty() ? -1 : mNucleotidesByLoci.lastKey(); }

    public NavigableMap<Integer, AminoAcid> aminoAcidsByLoci() { return mAminoAcidsByLoci; }
    public int minAminoAcidLocus() { return mAminoAcidsByLoci.isEmpty() ? -1 : mAminoAcidsByLoci.firstKey(); }
    public int maxAminoAcidLocus() { return mAminoAcidsByLoci.isEmpty() ? -1 : mAminoAcidsByLoci.lastKey(); }

    public void removeLowQualBases()
    {
        if(mIsQualFiltered)
            return;

        mIsQualFiltered = true;

        // cull any low-qual bases (they are retained in the raw bases)
        List<Integer> lociToRemove = Lists.newArrayList();
        for(Map.Entry<Integer, Nucleotide> entry : mNucleotidesByLoci.entrySet())
        {
            if(entry.getValue().isLowQual())
                lociToRemove.add(entry.getKey());
        }

        lociToRemove.forEach(mNucleotidesByLoci::remove);
    }

    public void buildAminoAcids()
    {
        ++mAminoAcidConversionCount;

        // build an amino-acid fragment from this nucleotide fragment, by creating amino acids for any complete codon
        mAminoAcidsByLoci.clear();

        for(Map.Entry<Integer, Nucleotide> entry : mNucleotidesByLoci.entrySet())
        {
            int locus = entry.getKey();
            if(!isCodonMultiple(locus))
                continue;

            if(mNucleotidesByLoci.containsKey(locus + 1) && mNucleotidesByLoci.containsKey(locus + 2))
            {
                int aaLocus = locus / 3;
                mAminoAcidsByLoci.put(aaLocus, new AminoAcid(aaLocus, formCodonAminoAcid(aaLocus)));
            }
        }
    }

    private String formCodonAminoAcid(int locus)
    {
        return FragmentUtils.formCodonAminoAcid(locus, mNucleotidesByLoci);
    }

    public boolean containsAminoAcidLoci(final Collection<Integer> loci)
    {
        return loci.stream().allMatch(this::containsAminoAcidLocus);
    }

    public boolean containsAminoAcidLocus(int locus)
    {
        return mAminoAcidsByLoci.containsKey(locus);
    }

    public String aminoAcid(int locus)
    {
        AminoAcid aminoAcid = mAminoAcidsByLoci.get(locus);
        if(aminoAcid == null)
            return "";

        return aminoAcid.acid();
    }

    public String aminoAcids(final Iterable<Integer> loci)
    {
        StringJoiner sj = new StringJoiner("");
        loci.forEach(x -> sj.add(aminoAcid(x)));
        return sj.toString();
    }

    public void filterAminoAcidsOnLoci(final Set<Integer> loci)
    {
        mAminoAcidsByLoci.keySet().removeIf(l -> !loci.contains(l));
    }

    public String getRawNucleotide(int locus)
    {
        Nucleotide nucleotide = mRawNucleotidesByLoci.get(locus);
        if(nucleotide == null)
            return "";

        return nucleotide.bases();
    }

    public String getLowQualAminoAcid(int locus)
    {
        int startNucleotideLocus = locus * 3;
        if(!mRawNucleotidesByLoci.containsKey(startNucleotideLocus))
            return "";

        if(!mRawNucleotidesByLoci.containsKey(startNucleotideLocus + 1))
            return "";

        if(!mRawNucleotidesByLoci.containsKey(startNucleotideLocus + 2))
            return "";

        return FragmentUtils.formCodonAminoAcid(locus, mRawNucleotidesByLoci);
    }

    public FragmentScope scope() { return mScope; }
    public void clearScope() { mScope = UNSET; }
    public boolean isScopeSet() { return mScope != UNSET; }

    public void setScope(final FragmentScope scope) { setScope(scope, false); }

    public void setScope(final FragmentScope scope, boolean override)
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

    @Override
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
    public void setAminoAcids(final Iterable<AminoAcid> aminoAcids)
    {
        for(AminoAcid aminoAcid : aminoAcids)
            mAminoAcidsByLoci.put(aminoAcid.locus(), aminoAcid);
    }
}
