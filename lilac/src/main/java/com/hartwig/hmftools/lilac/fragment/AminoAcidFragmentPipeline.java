package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
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
    private final List<AminoAcidFragment> mHighQualityTumorFragments;
    private final LilacConfig mConfig;
    private final List<NucleotideFragment> mReferenceFragments;

    public AminoAcidFragmentPipeline(
            final LilacConfig config, final List<NucleotideFragment> referenceFragments, final List<NucleotideFragment> tumorFragments)
    {
        mConfig = config;
        mReferenceFragments = referenceFragments;
        mMinBaseQuality = mConfig.MinBaseQual;
        mMinEvidence = mConfig.MinEvidence;
        mAminoAcidEnricher = new AminoAcidQualEnrichment(mMinEvidence);
        mNucleotideQualEnrichment = new NucleotideQualEnrichment(mMinBaseQuality, mMinEvidence);
        mHighQualityTumorFragments = qualFiltered(mMinBaseQuality, tumorFragments);
    }

    public final List<AminoAcidFragment> referencePhasingFragments(final HlaContext context)
    {
        String gene = "HLA-" + context.Gene;

        List<NucleotideFragment> geneReferenceFragments = mReferenceFragments.stream()
                .filter(x -> x.containsGene(gene)).collect(Collectors.toList());

        /*
        for(NucleotideFragment frag : geneReferenceFragments)
        {
            LL_LOGGER.info("read({}) has gene({})", frag.getId(), gene);
        }

        // NucleotideFragment frag1 = mReferenceFragments.stream().filter(x -> x.getId().equals("ST-E00285:220:HJG7JCCXY:6:1214:31132:40723")).findFirst().orElse(null);
         */

        List<AminoAcidFragment> referenceAminoAcids = process(context.AminoAcidBoundaries, geneReferenceFragments);

        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, nucFragments(referenceAminoAcids));
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, referenceAminoAcids);

        referenceAminoAcidCounts.writeVertically(mConfig.OutputFilePrefix + '.' + gene + ".aminoacids.txt");
        referenceNucleotideCounts.writeVertically(mConfig.OutputFilePrefix + '.' + gene + ".nucleotides.txt");

        return referenceAminoAcids;
    }

    public final List<AminoAcidFragment> referenceCoverageFragments()
    {
        return qualFiltered(mMinBaseQuality, mReferenceFragments);
    }

    public final List<AminoAcidFragment> tumorCoverageFragments()
    {
        if(mHighQualityTumorFragments.isEmpty())
            return Lists.newArrayList();

        List<AminoAcidFragment> referenceAminoAcids = referenceCoverageFragments();

        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, nucFragments(referenceAminoAcids));
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, referenceAminoAcids);

        SequenceCount tumorNucleotideCounts = SequenceCount.nucleotides(mMinEvidence, nucFragments(mHighQualityTumorFragments));
        SequenceCount tumorAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, mHighQualityTumorFragments);

        final List<SequenceCountDiff> nucleotideDifferences = SequenceCountDiff.create(referenceNucleotideCounts, tumorNucleotideCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<SequenceCountDiff> aminoAcidDifferences = SequenceCountDiff.create(referenceAminoAcidCounts, tumorAminoAcidCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<AminoAcidFragment> variantFilteredTumorAminoAcids = mHighQualityTumorFragments.stream()
                .filter(x -> !containsVariant(x, nucleotideDifferences, aminoAcidDifferences)).collect(
                Collectors.toList());

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

    private List<AminoAcidFragment> qualFiltered(int minBaseQuality, List<NucleotideFragment> fragments)
    {
        return fragments.stream()
                .map(x -> x.qualityFilter(minBaseQuality))
                .filter(x -> x.isNotEmpty())
                .map(x -> x.toAminoAcidFragment()).collect(
                Collectors.toList());
    }

    private List<AminoAcidFragment> process(final List<Integer> boundaries, final List<NucleotideFragment> fragments)
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
