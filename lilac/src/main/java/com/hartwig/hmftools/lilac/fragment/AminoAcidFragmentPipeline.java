package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.BASE_QUAL_FILTERED;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public class AminoAcidFragmentPipeline
{
    private final LilacConfig mConfig;

    private final int mMinEvidenceSupport_;
    private final double mMinEvidenceFactor_;
    private final double mMinHighQualEvidenceFactor_;

    private final List<Fragment> mHighQualRefAminoAcidFragments; // generated from input ref fragments, filtered and amino acids built

    // per-gene counts of bases and amino-acids with sufficient support
    private final Map<String, SequenceCount> mRefNucleotideCounts_;
    private final Map<String, SequenceCount> mRefAminoAcidCounts_;

    private final List<Fragment> mOriginalRefFragments; // copied and used for phasing only, will remain unchanged

    public AminoAcidFragmentPipeline(final LilacConfig config, final Collection<Fragment> referenceFragments)
    {
        mConfig = config;
        mOriginalRefFragments = referenceFragments.stream().map(FragmentUtils::copyNucleotideFragment).collect(Collectors.toList());

        mHighQualRefAminoAcidFragments = createHighQualAminoAcidFragments(referenceFragments);

        int fragmentCount = mHighQualRefAminoAcidFragments.size();

        mMinEvidenceSupport_ = config.MinEvidenceSupport_;
        mMinEvidenceFactor_ = config.MinEvidenceFactor_;
        mMinHighQualEvidenceFactor_ = config.MinHighQualEvidenceFactor_;

        mRefNucleotideCounts_ = Maps.newHashMap();
        mRefAminoAcidCounts_ = Maps.newHashMap();
    }

    public List<Fragment> highQualRefFragments() { return mHighQualRefAminoAcidFragments; }
    public int minEvidenceSupport_() { return mMinEvidenceSupport_; }
    public double minEvidenceFactor_() { return mMinEvidenceFactor_; }
    public double minHighQualEvidenceFactor_() { return mMinHighQualEvidenceFactor_; }

    public Map<String, SequenceCount> getReferenceAminoAcidCounts_()
    {
        // TODO:
        if(true)
        {
            System.out.println("");
        }

        // TODO:
        return Collections.unmodifiableMap(mRefAminoAcidCounts_);
    }

    private static List<Fragment> createHighQualAminoAcidFragments(final Iterable<Fragment> fragments)
    {
        List<Fragment> aminoAcidFragments = Lists.newArrayList();

        for(Fragment fragment : fragments)
        {
            fragment.removeLowQualBases();

            if(!fragment.hasNucleotides())
            {
                fragment.setScope(BASE_QUAL_FILTERED);
                continue;
            }

            fragment.buildAminoAcids();
            aminoAcidFragments.add(fragment);
        }

        return aminoAcidFragments;
    }

    public List<Set<String>> getReferenceNucleotides_()
    {
        // TODO:
        if(true)
        {
            System.out.println("");
        }

        return createSequenceSets_(mRefNucleotideCounts_);
    }

    private static List<Set<String>> createSequenceSets_(final Map<String, SequenceCount> countsMap)
    {
        NavigableMap<Integer, Set<String>> refAminoAcids = Maps.newTreeMap();
        for(SequenceCount seqCounts : countsMap.values())
        {
            for(int locus : seqCounts.seqCountsByLoci_().keySet())
            {
                refAminoAcids.computeIfAbsent(locus, l -> Sets.newHashSet());
                refAminoAcids.get(locus).addAll(seqCounts.seqCountsByLoci_().getOrDefault(locus, HashMultiset.create()).elementSet());
            }
        }

        return Lists.newArrayList(refAminoAcids.values());
    }

    public List<Fragment> referencePhasingFragments(final HlaContext context)
    {
        // filter and enrich fragments for phasing - all per gene

        // NOTE: this will create a copy of the fragment which will only be used in the scope of phasing and candidate selection,
        // not for anything to do with coverage
        String gene = context.geneName();

        // start with the unfiltered fragments again
        List<Fragment> geneRefNucFrags = mOriginalRefFragments.stream().filter(x -> x.containsGene(gene)).toList();
        List<Fragment> _geneRefNucFrags = geneRefNucFrags.stream().map(FragmentUtils::copyNucleotideFragment).toList();

        if(geneRefNucFrags.isEmpty())
            return Collections.emptyList();

        List<Fragment> highQualFrags = geneRefNucFrags.stream().map(FragmentUtils::copyNucleotideFragment).toList();
        highQualFrags.forEach(Fragment::removeLowQualBases);
        highQualFrags = highQualFrags.stream().filter(Fragment::hasNucleotides).toList();

        List<Fragment> qualEnrichedNucFrags = NucleotideFragmentQualEnrichment.qualityFilterFragments(
                mMinEvidenceSupport_, mMinEvidenceFactor_, mMinHighQualEvidenceFactor_, geneRefNucFrags, highQualFrags);

        int maxCommonAminoAcidExonBoundary = GENE_CACHE.MaxCommonAminoAcidExonBoundary;

        Set<Integer> aminoAcidBoundaries = context.AminoAcidBoundaries.stream()
                .filter(x -> x <= maxCommonAminoAcidExonBoundary).collect(Collectors.toSet());

        NucleotideSpliceEnrichment spliceEnricher = new NucleotideSpliceEnrichment(mMinEvidenceSupport_, mMinEvidenceFactor_, aminoAcidBoundaries);
        List<Fragment> spliceEnrichedNucFrags = spliceEnricher.applySpliceInfo(qualEnrichedNucFrags, highQualFrags);

        List<Fragment> enrichedAminoAcidFrags = AminoAcidQualEnrichment.qualityFilterAminoAcidFragments(mConfig, spliceEnrichedNucFrags, mMinEvidenceSupport_, mMinEvidenceFactor_);

        // cache support at each base and amino acid for later writing and recovery of low-qual support
        SequenceCount refNucleotideCounts = SequenceCount.nucleotides_(mMinEvidenceSupport_, mMinEvidenceFactor_, enrichedAminoAcidFrags);
        SequenceCount refAminoAcidCounts = SequenceCount.aminoAcids_(mMinEvidenceSupport_, mMinEvidenceFactor_, enrichedAminoAcidFrags);

        // TODO:
        if(gene.equals("HLA-B"))
        {
            List<Fragment> _frags = mOriginalRefFragments.stream().filter(x -> x.readGene().equals("HLA-B")).sorted(Comparator.comparingInt((Fragment x) -> x.reads().get(0).bamRecord().getAlignmentStart()).reversed()).toList();

            SequenceCount _refNucleotideCounts = SequenceCount.nucleotides_(mMinEvidenceSupport_, mMinEvidenceFactor_, _geneRefNucFrags);
            SequenceCount _fragsNucCounts = SequenceCount.nucleotides_(mMinEvidenceSupport_, mMinEvidenceFactor_, _frags);

            _geneRefNucFrags.forEach(Fragment::buildAminoAcids);
            _frags.forEach(Fragment::buildAminoAcids);

            SequenceCount _refAminoAcidCounts = SequenceCount.aminoAcids_(mMinEvidenceSupport_, mMinEvidenceFactor_, _geneRefNucFrags);
            SequenceCount _fragsAminoCounts = SequenceCount.aminoAcids_(mMinEvidenceSupport_, mMinEvidenceFactor_, _frags);

            System.out.println("");
        }

        setCounts_(gene, refNucleotideCounts, refAminoAcidCounts);

        return enrichedAminoAcidFrags;
    }

    private synchronized void setCounts_(final String gene, final SequenceCount refNucleotideCounts, final SequenceCount refAminoAcidCounts)
    {
        // TODO: we are setting thins here
        mRefNucleotideCounts_.put(gene, refNucleotideCounts);
        mRefAminoAcidCounts_.put(gene, refAminoAcidCounts);
    }

    public List<Fragment> calcComparisonCoverageFragments(final Iterable<Fragment> comparisonFragments)
    {
        List<Fragment> highQualFragments = createHighQualAminoAcidFragments(comparisonFragments);

        if(highQualFragments.isEmpty())
            return Lists.newArrayList();

        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides_(mMinEvidenceSupport_, mMinEvidenceFactor_, mHighQualRefAminoAcidFragments);
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids_(mMinEvidenceSupport_, mMinEvidenceFactor_, mHighQualRefAminoAcidFragments);

        SequenceCount tumorNucleotideCounts = SequenceCount.nucleotides_(mMinEvidenceSupport_, mMinEvidenceFactor_, highQualFragments);
        SequenceCount tumorAminoAcidCounts = SequenceCount.aminoAcids_(mMinEvidenceSupport_, mMinEvidenceFactor_, highQualFragments);

        final List<SequenceCountDiff> nucleotideDifferences = SequenceCountDiff.create(referenceNucleotideCounts, tumorNucleotideCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<SequenceCountDiff> aminoAcidDifferences = SequenceCountDiff.create(referenceAminoAcidCounts, tumorAminoAcidCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<Fragment> variantFilteredTumorAminoAcids = highQualFragments.stream()
                .filter(x -> !containsVariant(x, nucleotideDifferences, aminoAcidDifferences)).collect(Collectors.toList());

        return variantFilteredTumorAminoAcids;
    }

    private static boolean containsVariant(final Fragment fragment, final Collection<SequenceCountDiff> nucelotideVariants,
            final Collection<SequenceCountDiff> aminoAcidVariants)
    {
        return nucelotideVariants.stream().anyMatch(x -> containsNucleotideVariant(fragment, x))
                || aminoAcidVariants.stream().anyMatch(x -> containsAminoAcidVariant(fragment, x));
    }

    private static boolean containsNucleotideVariant(final Fragment fragment, final SequenceCountDiff variant)
    {
        return fragment.containsNucleotideLocus(variant.Loci) && fragment.nucleotide(variant.Loci).equals(variant.Sequence);
    }

    private static boolean containsAminoAcidVariant(final Fragment fragment, final SequenceCountDiff variant)
    {
        return fragment.containsAminoAcidLocus(variant.Loci) && fragment.aminoAcid(variant.Loci).equals(variant.Sequence);
    }

    public void writeCounts(final LilacConfig config)
    {
        for(Map.Entry<String, SequenceCount> entry : mRefAminoAcidCounts_.entrySet())
        {
            String gene = entry.getKey();
            entry.getValue().writeVertically(config.formFileId(gene + ".aminoacids.txt"), gene);
        }

        for(Map.Entry<String, SequenceCount> entry : mRefNucleotideCounts_.entrySet())
        {
            String gene = entry.getKey();
            entry.getValue().writeVertically(config.formFileId(gene + ".nucleotides.txt"), gene);
        }
    }
}
