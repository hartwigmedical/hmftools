package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.MAX_AMINO_ACID_BOUNDARY;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragment.nucFragments;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.hartwig.hmftools.lilac.SequenceCountDiff;
import com.hartwig.hmftools.lilac.hla.HlaContext;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class AminoAcidFragmentPipeline
{
    private final int mMinBaseQuality;
    private final int mMinEvidence;

    private final List<NucleotideFragment> mRefNucleotideFragments;

    private final List<AminoAcidFragment> mHighQualRefAminoAcidFragments;
    private final List<AminoAcidFragment> mHighQualTumorFragments;

    private final Map<String,SequenceCount> mRefNucleotideCounts;
    private final Map<String,SequenceCount> mRefAminoAcidCounts;

    private final AminoAcidQualEnrichment mAminoAcidEnricher;
    private final NucleotideQualEnrichment mNucleotideQualEnrichment;

    public AminoAcidFragmentPipeline(
            final LilacConfig config, final List<NucleotideFragment> referenceFragments, final List<NucleotideFragment> tumorFragments)
    {
        mRefNucleotideFragments = referenceFragments;
        mMinBaseQuality = config.MinBaseQual;
        mMinEvidence = config.MinEvidence;

        mHighQualRefAminoAcidFragments = createHighQualAminoAcidFragments(mRefNucleotideFragments);
        mHighQualTumorFragments = createHighQualAminoAcidFragments(tumorFragments);

        // only used for phasing fragments
        mAminoAcidEnricher = new AminoAcidQualEnrichment(mMinEvidence);
        mNucleotideQualEnrichment = new NucleotideQualEnrichment(mMinBaseQuality, mMinEvidence);

        mRefNucleotideCounts = Maps.newHashMap();
        mRefAminoAcidCounts = Maps.newHashMap();
    }

    public List<AminoAcidFragment> referenceAminoAcidFragments() { return mHighQualRefAminoAcidFragments; }

    // public Map<String,SequenceCount> referenceNucleotideCounts() { return mRefNucleotideCounts; }
    public Map<String,SequenceCount> referenceAminoAcidCounts() { return mRefAminoAcidCounts; }

    private List<AminoAcidFragment> createHighQualAminoAcidFragments(final List<NucleotideFragment> fragments)
    {
        List<AminoAcidFragment> aminoAcidFragments = Lists.newArrayList();

        for(NucleotideFragment fragment : fragments)
        {
            NucleotideFragment highQualFragment = fragment.qualityFilter(mMinBaseQuality);

            if(highQualFragment.isEmpty())
                continue;

            AminoAcidFragment aaFragment = highQualFragment.toAminoAcidFragment();
            aaFragment.setAllQualitytNucleotideFragment(fragment);
            aminoAcidFragments.add(aaFragment);
        }

        return aminoAcidFragments;
    }

    public List<AminoAcidFragment> referencePhasingFragments(final HlaContext context)
    {
        String gene = context.geneName();

        // start with the unfiltered fragments again
        List<NucleotideFragment> geneReferenceFragments = mRefNucleotideFragments.stream()
                .filter(x -> x.containsGene(gene)).collect(Collectors.toList());

        List<AminoAcidFragment> refAminoAcids = applyQualAndSpliceChecks(context.AminoAcidBoundaries, geneReferenceFragments);

        // cache support at each base and amino acid for later writing and recovery of low-qual support
        SequenceCount refNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, nucFragments(refAminoAcids));
        mRefNucleotideCounts.put(gene, refNucleotideCounts);

        SequenceCount refAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, refAminoAcids);
        mRefAminoAcidCounts.put(gene, refAminoAcidCounts);

        return refAminoAcids;
    }

    public List<Set<String>> getReferenceAminoAcids()
    {
        List<Set<String>> refAminoAcids = Lists.newArrayList();

        for(SequenceCount seqCounts : mRefAminoAcidCounts.values())
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

    private List<AminoAcidFragment> applyQualAndSpliceChecks(final List<Integer> boundaries, final List<NucleotideFragment> fragments)
    {
        if(fragments.isEmpty())
            return Lists.newArrayList();

        List<NucleotideFragment> qualEnriched = mNucleotideQualEnrichment.enrich(fragments);

        NucleotideSpliceEnrichment spliceEnricher = new NucleotideSpliceEnrichment(
                mMinBaseQuality, mMinEvidence, boundaries.stream().filter(x -> x <= MAX_AMINO_ACID_BOUNDARY).collect(Collectors.toSet()));

        List<NucleotideFragment> spliceEnriched = spliceEnricher.enrich(qualEnriched);

        return mAminoAcidEnricher.enrich(spliceEnriched);
    }

    public List<AminoAcidFragment> tumorCoverageFragments()
    {
        if(mHighQualTumorFragments.isEmpty())
            return Lists.newArrayList();

        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, nucFragments(mHighQualRefAminoAcidFragments));
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, mHighQualRefAminoAcidFragments);

        SequenceCount tumorNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, nucFragments(mHighQualTumorFragments));
        SequenceCount tumorAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, mHighQualTumorFragments);

        final List<SequenceCountDiff> nucleotideDifferences = SequenceCountDiff.create(referenceNucleotideCounts, tumorNucleotideCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<SequenceCountDiff> aminoAcidDifferences = SequenceCountDiff.create(referenceAminoAcidCounts, tumorAminoAcidCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<AminoAcidFragment> variantFilteredTumorAminoAcids = mHighQualTumorFragments.stream()
                .filter(x -> !containsVariant(x, nucleotideDifferences, aminoAcidDifferences)).collect(Collectors.toList());

        return variantFilteredTumorAminoAcids;
    }

    private static boolean containsVariant(
            final AminoAcidFragment fragment,
            final List<SequenceCountDiff> nucelotideVariants, final List<SequenceCountDiff> aminoAcidVariants)
    {
        return nucelotideVariants.stream().anyMatch(x -> containsNucleotideVariant(fragment, x))
                || aminoAcidVariants.stream().anyMatch(x -> containsAminoAcidVariant(fragment, x));
    }

    private static boolean containsNucleotideVariant(final AminoAcidFragment fragment, SequenceCountDiff variant)
    {
        return fragment.containsNucleotide(variant.Loci) && fragment.nucleotide(variant.Loci).equals(variant.Sequence);
    }

    private static boolean containsAminoAcidVariant(final AminoAcidFragment fragment, SequenceCountDiff variant)
    {
        return fragment.containsAminoAcid(variant.Loci) && fragment.aminoAcid(variant.Loci).equals(variant.Sequence);
    }

}
