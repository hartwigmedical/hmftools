package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.BASE_QUAL_FILTERED;

import java.util.Collection;
import java.util.Collections;
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

    private final double mMinEvidenceFactor;
    private final double mMinHighQualEvidenceFactor;
    private final int mMinVafFilterDepth;

    private final List<Fragment> mHighQualRefAminoAcidFragments; // generated from input ref fragments, filtered and amino acids built

    // per-gene counts of bases and amino-acids with sufficient support
    private final Map<String, SequenceCount> mRefNucleotideCounts;
    private final Map<String, SequenceCount> mRefAminoAcidCounts;

    private final List<Fragment> mOriginalRefFragments; // copied and used for phasing only, will remain unchanged

    public AminoAcidFragmentPipeline(final LilacConfig config, final Collection<Fragment> referenceFragments)
    {
        mConfig = config;
        mOriginalRefFragments = referenceFragments.stream().map(FragmentUtils::copyNucleotideFragment).collect(Collectors.toList());

        mHighQualRefAminoAcidFragments = createHighQualAminoAcidFragments(referenceFragments);

        mMinEvidenceFactor = config.MinEvidenceFactor;
        mMinHighQualEvidenceFactor = config.MinHighQualEvidenceFactor;
        mMinVafFilterDepth = config.MinVafFilterDepth;

        mRefNucleotideCounts = Maps.newHashMap();
        mRefAminoAcidCounts = Maps.newHashMap();
    }

    public List<Fragment> highQualRefFragments()
    {
        return mHighQualRefAminoAcidFragments;
    }

    public Map<String, SequenceCount> getReferenceAminoAcidCounts()
    {
        return mRefAminoAcidCounts;
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

    public List<Set<String>> getReferenceNucleotides()
    {
        return createSequenceSets(mRefNucleotideCounts);
    }

    private static List<Set<String>> createSequenceSets(final Map<String, SequenceCount> countsMap)
    {
        NavigableMap<Integer, Set<String>> refAminoAcids = Maps.newTreeMap();
        for(SequenceCount seqCounts : countsMap.values())
        {
            for(int locus : seqCounts.seqCountsByLoci().keySet())
            {
                refAminoAcids.computeIfAbsent(locus, l -> Sets.newHashSet());
                refAminoAcids.get(locus).addAll(seqCounts.seqCountsByLoci().getOrDefault(locus, HashMultiset.create()).elementSet());
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
        if(geneRefNucFrags.isEmpty())
        {
            return Collections.emptyList();
        }

        List<Fragment> highQualFrags = geneRefNucFrags.stream().map(FragmentUtils::copyNucleotideFragment).toList();
        highQualFrags.forEach(Fragment::removeLowQualBases);
        highQualFrags = highQualFrags.stream().filter(Fragment::hasNucleotides).toList();

        List<Fragment> qualEnrichedNucFrags = NucleotideFragmentQualEnrichment.qualityFilterFragments(
                mMinVafFilterDepth, mMinEvidenceFactor, mMinHighQualEvidenceFactor, geneRefNucFrags, highQualFrags);

        int maxCommonAminoAcidExonBoundary = GENE_CACHE.MaxCommonAminoAcidExonBoundary;

        Set<Integer> aminoAcidBoundaries = context.AminoAcidBoundaries.stream()
                .filter(x -> x <= maxCommonAminoAcidExonBoundary).collect(Collectors.toSet());

        NucleotideSpliceEnrichment spliceEnricher = new NucleotideSpliceEnrichment(
                mMinVafFilterDepth, mMinEvidenceFactor, aminoAcidBoundaries);
        List<Fragment> spliceEnrichedNucFrags = spliceEnricher.applySpliceInfo(qualEnrichedNucFrags, highQualFrags);

        List<Fragment> enrichedAminoAcidFrags = AminoAcidQualEnrichment.qualityFilterAminoAcidFragments(mConfig, spliceEnrichedNucFrags);

        // cache support at each base and amino acid for later writing and recovery of low-qual support
        SequenceCount refNucleotideCounts = SequenceCount.nucleotides(mMinVafFilterDepth, mMinEvidenceFactor, enrichedAminoAcidFrags);
        SequenceCount refAminoAcidCounts = SequenceCount.aminoAcids(mMinVafFilterDepth, mMinEvidenceFactor, enrichedAminoAcidFrags);
        setCounts(gene, refNucleotideCounts, refAminoAcidCounts);

        return enrichedAminoAcidFrags;
    }

    private synchronized void setCounts(final String gene, final SequenceCount refNucleotideCounts, final SequenceCount refAminoAcidCounts)
    {
        mRefNucleotideCounts.put(gene, refNucleotideCounts);
        mRefAminoAcidCounts.put(gene, refAminoAcidCounts);
    }

    public List<Fragment> calcComparisonCoverageFragments(final Iterable<Fragment> comparisonFragments)
    {
        List<Fragment> highQualFragments = createHighQualAminoAcidFragments(comparisonFragments);

        if(highQualFragments.isEmpty())
        {
            return Lists.newArrayList();
        }

        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides(
                mMinVafFilterDepth, mMinEvidenceFactor, mHighQualRefAminoAcidFragments);
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(
                mMinVafFilterDepth, mMinEvidenceFactor, mHighQualRefAminoAcidFragments);

        SequenceCount tumorNucleotideCounts = SequenceCount.nucleotides(mMinVafFilterDepth, mMinEvidenceFactor, highQualFragments);
        SequenceCount tumorAminoAcidCounts = SequenceCount.aminoAcids(mMinVafFilterDepth, mMinEvidenceFactor, highQualFragments);

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
        for(Map.Entry<String, SequenceCount> entry : mRefAminoAcidCounts.entrySet())
        {
            String gene = entry.getKey();
            entry.getValue().writeVertically(config.formFileId(gene + ".aminoacids.txt"));
        }

        for(Map.Entry<String, SequenceCount> entry : mRefNucleotideCounts.entrySet())
        {
            String gene = entry.getKey();
            entry.getValue().writeVertically(config.formFileId(gene + ".nucleotides.txt"));
        }
    }
}
