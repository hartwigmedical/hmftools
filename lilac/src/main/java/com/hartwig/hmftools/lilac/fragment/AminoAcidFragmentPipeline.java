package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.MAX_AMINO_ACID_BOUNDARY;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragment.nucFragments;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.hartwig.hmftools.lilac.SequenceCountDiff;
import com.hartwig.hmftools.lilac.hla.HlaContext;

import java.util.List;
import java.util.stream.Collectors;

public class AminoAcidFragmentPipeline
{
    private final int mMinBaseQuality;
    private final int mMinEvidence;
    private final AminoAcidQualEnrichment mAminoAcidEnricher;
    private final NucleotideQualEnrichment mNucleotideQualEnrichment;
    private final List<AminoAcidFragment> mHighQualTumorFragments;
    private final LilacConfig mConfig;

    private final List<NucleotideFragment> mReferenceFragments;
    private final List<AminoAcidFragment> mHighQualReferenceFragments;

    public AminoAcidFragmentPipeline(
            final LilacConfig config, final List<NucleotideFragment> referenceFragments, final List<NucleotideFragment> tumorFragments)
    {
        mConfig = config;
        mReferenceFragments = referenceFragments;
        mMinBaseQuality = mConfig.MinBaseQual;
        mMinEvidence = mConfig.MinEvidence;

        mAminoAcidEnricher = new AminoAcidQualEnrichment(mMinEvidence);
        mNucleotideQualEnrichment = new NucleotideQualEnrichment(mMinBaseQuality, mMinEvidence);

        mHighQualReferenceFragments = applyQualityFilter(mReferenceFragments);
        mHighQualTumorFragments = applyQualityFilter(tumorFragments);
    }

    public final List<AminoAcidFragment> referencePhasingFragments(final HlaContext context)
    {
        String gene = context.geneName();

        // start with the unfiltered fragments again
        List<NucleotideFragment> geneReferenceFragments = mReferenceFragments.stream()
                .filter(x -> x.containsGene(gene)).collect(Collectors.toList());

        List<AminoAcidFragment> referenceAminoAcids = applyQualAndSpliceChecks(context.AminoAcidBoundaries, geneReferenceFragments);

        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, nucFragments(referenceAminoAcids));
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, referenceAminoAcids);

        referenceAminoAcidCounts.writeVertically(mConfig.outputPrefix() + '.' + gene + ".aminoacids.txt");
        referenceNucleotideCounts.writeVertically(mConfig.outputPrefix() + '.' + gene + ".nucleotides.txt");

        return referenceAminoAcids;
    }

    public final List<AminoAcidFragment> referenceCoverageFragments() { return mHighQualReferenceFragments; }

    public final List<AminoAcidFragment> tumorCoverageFragments()
    {
        if(mHighQualTumorFragments.isEmpty())
            return Lists.newArrayList();

        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, nucFragments(mHighQualReferenceFragments));
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, mHighQualReferenceFragments);

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

    private List<AminoAcidFragment> applyQualityFilter(List<NucleotideFragment> fragments)
    {
        return fragments.stream()
                .map(x -> x.qualityFilter(mMinBaseQuality))
                .filter(x -> x.isNotEmpty())
                .map(x -> x.toAminoAcidFragment()).collect(
                Collectors.toList());
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
}
