package com.hartwig.hmftools.lilac.fragment;

import static java.lang.Math.max;

import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.BASE_QUAL_FILTERED;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.copyNucleotideFragment;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public class AminoAcidFragmentPipeline
{
    private final LilacConfig mConfig;
    private final double mMinEvidence;
    private final double mMinHighQualEvidence;

    private final List<Fragment> mHighQualRefAminoAcidFragments; // generated from input ref fragments, filtered and amino acids built

    // per-gene counts of bases and amino-acids with sufficient support
    private final Map<String,SequenceCount> mRefNucleotideCounts;
    private final Map<String,SequenceCount> mRefAminoAcidCounts;

    private final List<Fragment> mOriginalRefFragments; // copied and used for phasing only, will remain unchanged

    public AminoAcidFragmentPipeline(final LilacConfig config, final List<Fragment> referenceFragments)
    {
        mConfig = config;
        mOriginalRefFragments = referenceFragments.stream().map(x -> copyNucleotideFragment(x)).collect(Collectors.toList());

        mHighQualRefAminoAcidFragments = createHighQualAminoAcidFragments(referenceFragments);

        int fragmentCount = mHighQualRefAminoAcidFragments.size();

        mMinEvidence = max(config.MinEvidence, fragmentCount * config.MinNucleotideEvidenceFactor);
        mMinHighQualEvidence = fragmentCount * config.MinNucleotideHighQualEvidenceFactor;

        mRefNucleotideCounts = Maps.newHashMap();
        mRefAminoAcidCounts = Maps.newHashMap();
    }

    public List<Fragment> highQualRefFragments() { return mHighQualRefAminoAcidFragments; }
    public double minEvidence() { return mMinEvidence; }
    public double minHighQualEvidence() { return mMinHighQualEvidence; }

    public Map<String,SequenceCount> getReferenceAminoAcidCounts() { return mRefAminoAcidCounts; }

    private List<Fragment> createHighQualAminoAcidFragments(final List<Fragment> fragments)
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

    private static List<Set<String>> createSequenceSets(final Map<String,SequenceCount> countsMap)
    {
        List<Set<String>> refAminoAcids = Lists.newArrayList();

        for(SequenceCount seqCounts : countsMap.values())
        {
            for(int locus = 0; locus < seqCounts.getLength(); ++locus)
            {
                Set<String> aminoAcids;
                if(locus >= refAminoAcids.size())
                {
                    aminoAcids = Sets.newHashSet();
                    refAminoAcids.add(locus, aminoAcids);
                }
                else
                {
                    aminoAcids = refAminoAcids.get(locus);
                }

                aminoAcids.addAll(seqCounts.get(locus).keySet());
            }
        }

        return refAminoAcids;
    }

    public List<Fragment> referencePhasingFragments(final HlaContext context)
    {
        // filter and enrich fragments for phasing - all per gene

        // NOTE: this will create a copy of the fragment which will only be used in the scope of phasing and candidate selection,
        // not for anything to do with coverage
        String gene = context.geneName();

        // start with the unfiltered fragments again
        List<Fragment> geneRefNucFrags = mOriginalRefFragments.stream()
                .filter(x -> x.containsGene(gene))
                .collect(Collectors.toList());

        if(geneRefNucFrags.isEmpty())
            return Collections.emptyList();

        List<Fragment> highQualFrags = geneRefNucFrags.stream()
                .map(x -> copyNucleotideFragment(x))
                .collect(Collectors.toList());

        highQualFrags.forEach(x -> x.removeLowQualBases());
        highQualFrags = highQualFrags.stream().filter(x -> x.hasNucleotides()).collect(Collectors.toList());

        List<Fragment> qualEnrichedNucFrags = NucleotideFragmentQualEnrichment.qualityFilterFragments(
                mMinEvidence, mMinHighQualEvidence, geneRefNucFrags, highQualFrags);

        int maxCommonAminoAcidExonBoundary = GENE_CACHE.MaxCommonAminoAcidExonBoundary;

        Set<Integer> aminoAcidBoundaries = context.AminoAcidBoundaries.stream()
                .filter(x -> x <= maxCommonAminoAcidExonBoundary).collect(Collectors.toSet());

        NucleotideSpliceEnrichment spliceEnricher = new NucleotideSpliceEnrichment(mMinEvidence, aminoAcidBoundaries);
        List<Fragment> spliceEnrichedNucFrags = spliceEnricher.applySpliceInfo(qualEnrichedNucFrags, highQualFrags);

        List<Fragment> enrichedAminoAcidFrags = AminoAcidQualEnrichment.qualityFilterAminoAcidFragments(mConfig, spliceEnrichedNucFrags, mMinEvidence);

        // cache support at each base and amino acid for later writing and recovery of low-qual support
        SequenceCount refNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, enrichedAminoAcidFrags);
        SequenceCount refAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, enrichedAminoAcidFrags);
        setCounts(gene, refNucleotideCounts, refAminoAcidCounts);

        return enrichedAminoAcidFrags;
    }

    private synchronized void setCounts(final String gene, final SequenceCount refNucleotideCounts, final SequenceCount refAminoAcidCounts)
    {
        mRefNucleotideCounts.put(gene, refNucleotideCounts);
        mRefAminoAcidCounts.put(gene, refAminoAcidCounts);
    }

    public List<Fragment> calcComparisonCoverageFragments(final List<Fragment> comparisonFragments)
    {
        List<Fragment> highQualFragments = createHighQualAminoAcidFragments(comparisonFragments);

        if(highQualFragments.isEmpty())
            return Lists.newArrayList();

        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, mHighQualRefAminoAcidFragments);
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, mHighQualRefAminoAcidFragments);

        SequenceCount tumorNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, highQualFragments);
        SequenceCount tumorAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, highQualFragments);

        final List<SequenceCountDiff> nucleotideDifferences = SequenceCountDiff.create(referenceNucleotideCounts, tumorNucleotideCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<SequenceCountDiff> aminoAcidDifferences = SequenceCountDiff.create(referenceAminoAcidCounts, tumorAminoAcidCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<Fragment> variantFilteredTumorAminoAcids = highQualFragments.stream()
                .filter(x -> !containsVariant(x, nucleotideDifferences, aminoAcidDifferences)).collect(Collectors.toList());

        return variantFilteredTumorAminoAcids;
    }

    private static boolean containsVariant(
            final Fragment fragment, final List<SequenceCountDiff> nucelotideVariants, final List<SequenceCountDiff> aminoAcidVariants)
    {
        return nucelotideVariants.stream().anyMatch(x -> containsNucleotideVariant(fragment, x))
                || aminoAcidVariants.stream().anyMatch(x -> containsAminoAcidVariant(fragment, x));
    }

    private static boolean containsNucleotideVariant(final Fragment fragment, SequenceCountDiff variant)
    {
        return fragment.containsNucleotideLocus(variant.Loci) && fragment.nucleotide(variant.Loci).equals(variant.Sequence);
    }

    private static boolean containsAminoAcidVariant(final Fragment fragment, SequenceCountDiff variant)
    {
        return fragment.containsAminoAcidLocus(variant.Loci) && fragment.aminoAcid(variant.Loci).equals(variant.Sequence);
    }

    public void writeCounts(final LilacConfig config)
    {
        for(Map.Entry<String,SequenceCount> entry : mRefAminoAcidCounts.entrySet())
        {
            String gene = entry.getKey();
            entry.getValue().writeVertically(config.formFileId(gene + ".aminoacids.txt"));
        }

        for(Map.Entry<String,SequenceCount> entry : mRefNucleotideCounts.entrySet())
        {
            String gene = entry.getKey();
            entry.getValue().writeVertically(config.formFileId(gene + ".nucleotides.txt"));
        }
    }
}
