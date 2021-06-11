package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.MAX_AMINO_ACID_BOUNDARY;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.BASE_QUAL_FILTERED;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.copyNucleotideFragment;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.hla.HlaContext;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class AminoAcidFragmentPipeline
{
    private final int mMinBaseQuality;
    private final double mMinEvidence;
    private final double mMinHighQualEvidence;

    private final List<Fragment> mHighQualRefAminoAcidFragments;

    // per-gene counts of bases and amino-acids with sufficient support
    private final Map<String,SequenceCount> mRefNucleotideCounts;
    private final Map<String,SequenceCount> mRefAminoAcidCounts;

    private final List<Fragment> mRefNucPhasingFragments; // copied and used for phasing only

    public AminoAcidFragmentPipeline(final LilacConfig config, final List<Fragment> referenceFragments)
    {
        mMinBaseQuality = config.MinBaseQual;

        mRefNucPhasingFragments = referenceFragments.stream().map(x -> copyNucleotideFragment(x)).collect(Collectors.toList());

        mHighQualRefAminoAcidFragments = createHighQualAminoAcidFragments(referenceFragments);

        int fragmentCount = mHighQualRefAminoAcidFragments.size();

        mMinEvidence = config.calcMinEvidence(fragmentCount);
        mMinHighQualEvidence = config.calcMinHighQualEvidence(fragmentCount);

        mRefNucleotideCounts = Maps.newHashMap();
        mRefAminoAcidCounts = Maps.newHashMap();
    }

    public List<Fragment> getReferenceFragments() { return mHighQualRefAminoAcidFragments; }

    public Map<String,SequenceCount> getReferenceAminoAcidCounts() { return mRefAminoAcidCounts; }

    private List<Fragment> createHighQualAminoAcidFragments(final List<Fragment> fragments)
    {
        List<Fragment> aminoAcidFragments = Lists.newArrayList();

        for(Fragment fragment : fragments)
        {
            fragment.qualityFilter(mMinBaseQuality);

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
        List<Fragment> geneRefNucFrags = mRefNucPhasingFragments.stream()
                .filter(x -> x.containsGene(gene))
                .collect(Collectors.toList());

        if(geneRefNucFrags.isEmpty())
            return Lists.newArrayList();

        List<Fragment> highQualFrags = geneRefNucFrags.stream()
                .map(x -> copyNucleotideFragment(x))
                .collect(Collectors.toList());

        highQualFrags.forEach(x -> x.qualityFilter(mMinBaseQuality));
        highQualFrags = highQualFrags.stream().filter(x -> x.hasNucleotides()).collect(Collectors.toList());

        List<Fragment> qualEnrichedNucFrags = NucleotideFragmentQualEnrichment.enrich(
                mMinEvidence, mMinHighQualEvidence, geneRefNucFrags, highQualFrags);

        Set<Integer> aminoAcidBoundaries = context.AminoAcidBoundaries.stream()
                .filter(x -> x <= MAX_AMINO_ACID_BOUNDARY).collect(Collectors.toSet());

        NucleotideSpliceEnrichment spliceEnricher = new NucleotideSpliceEnrichment(mMinBaseQuality, mMinEvidence, aminoAcidBoundaries);
        List<Fragment> spliceEnrichedNucFrags = spliceEnricher.enrich(qualEnrichedNucFrags, highQualFrags);

        AminoAcidQualEnrichment aminoAcidEnricher = new AminoAcidQualEnrichment(mMinEvidence);
        List<Fragment> enrichedAminoAcidFrags = aminoAcidEnricher.enrich(spliceEnrichedNucFrags);

        // cache support at each base and amino acid for later writing and recovery of low-qual support
        SequenceCount refNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, enrichedAminoAcidFrags);
        mRefNucleotideCounts.put(gene, refNucleotideCounts);

        SequenceCount refAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, enrichedAminoAcidFrags);
        mRefAminoAcidCounts.put(gene, refAminoAcidCounts);

        return enrichedAminoAcidFrags;
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
        return fragment.containsNucleotide(variant.Loci) && fragment.nucleotide(variant.Loci).equals(variant.Sequence);
    }

    private static boolean containsAminoAcidVariant(final Fragment fragment, SequenceCountDiff variant)
    {
        return fragment.containsAminoAcid(variant.Loci) && fragment.aminoAcid(variant.Loci).equals(variant.Sequence);
    }

    public void writeCounts(final String outputPrefix)
    {
        for(Map.Entry<String,SequenceCount> entry : mRefAminoAcidCounts.entrySet())
        {
            String gene = entry.getKey();
            entry.getValue().writeVertically(outputPrefix + '.' + gene + ".aminoacids.txt");
        }

        for(Map.Entry<String,SequenceCount> entry : mRefNucleotideCounts.entrySet())
        {
            String gene = entry.getKey();
            entry.getValue().writeVertically(outputPrefix + '.' + gene + ".nucleotides.txt");
        }
    }

}
