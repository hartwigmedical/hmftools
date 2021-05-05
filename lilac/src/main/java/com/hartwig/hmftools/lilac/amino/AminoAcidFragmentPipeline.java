package com.hartwig.hmftools.lilac.amino;

import static com.hartwig.hmftools.lilac.LilacConstants.MAX_AMINO_ACID_BOUNDARY;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.hartwig.hmftools.lilac.SequenceCountDiff;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment;
import com.hartwig.hmftools.lilac.nuc.NucleotideQualEnrichment;
import com.hartwig.hmftools.lilac.nuc.NucleotideSpliceEnrichment;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
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
                .filter(x -> x.getGenes().contains(gene)).collect(Collectors.toList());

        List<AminoAcidFragment> referenceAminoAcids = process(context.AminoAcidBoundaries, geneReferenceFragments);

        // TODO, CHECK - need to cast to super type?
        SequenceCount referenceNucleotideCounts = null; // SequenceCount.nucleotides(mMinEvidence, referenceAminoAcids);
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

        // TODO: auto-casting from nuc frag to AA frag
        List<AminoAcidFragment> referenceAminoAcids = referenceCoverageFragments();
        SequenceCount referenceNucleotideCounts = null; // SequenceCount.nucleotides(mMinEvidence, referenceAminoAcids);
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, referenceAminoAcids);
        SequenceCount tumorNucleotideCounts = null; // SequenceCount.nucleotides(mMinEvidence, mHighQualityTumorFragments);
        SequenceCount tumorAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, mHighQualityTumorFragments);

        /*
                val nucleotideDifferences = SequenceCountDiff.create(referenceNucleotideCounts, tumorNucleotideCounts).filter { it.tumorCount > 0 }
        val aminoAcidDifferences = SequenceCountDiff.create(referenceAminoAcidCounts, tumorAminoAcidCounts).filter { it.tumorCount > 0 }

        val variantFilteredTumorAminoAcids = highQualityTumorFragments.filter { !it.containsVariant(nucleotideDifferences, aminoAcidDifferences) }
        return variantFilteredTumorAminoAcids

         */

        final List<SequenceCountDiff> nucleotideDifferences = SequenceCountDiff.create(referenceNucleotideCounts, tumorNucleotideCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<SequenceCountDiff> aminoAcidDifferences = SequenceCountDiff.create(referenceAminoAcidCounts, tumorAminoAcidCounts)
                .stream().filter(x -> x.TumorCount > 0).collect(Collectors.toList());

        final List<AminoAcidFragment> variantFilteredTumorAminoAcids = mHighQualityTumorFragments.stream()
                .filter(x -> !containsVariant(nucleotideDifferences, aminoAcidDifferences)).collect(
                Collectors.toList());

        return variantFilteredTumorAminoAcids;
    }

    private final boolean containsVariant(
            final List<SequenceCountDiff> nucelotideVariants, final List<SequenceCountDiff> aminoAcidVariants)
    {
        /*

        return nucelotideVariants.any { this.containsNucleotideVariant(it) } || aminoAcidVariants.any { this.containsAminoAcidVariant(it) }

         */
        return true;
    }

    private final boolean containsNucleotideVariant(final AminoAcidFragment fragment, SequenceCountDiff variant)
    {
        return fragment.containsNucleotide(variant.Loci) && fragment.nucleotide(variant.Loci).equals(variant.Sequence);
    }

    private final boolean containsAminoAcidVariant(final AminoAcidFragment fragment, SequenceCountDiff variant)
    {
        return fragment.containsAminoAcid(variant.Loci) && fragment.aminoAcid(variant.Loci).equals(variant.Sequence);
    }

    private final List<AminoAcidFragment> qualFiltered(int minBaseQuality, List<NucleotideFragment> fragments)
    {
        /*
        val qualityFilteredFragments = fragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }
        return qualityFilteredFragments.map { it.toAminoAcidFragment() }

         */

        return fragments.stream()
                .map(x -> x.qualityFilter(minBaseQuality))
                .filter(x -> x.isNotEmpty())
                .map(x -> x.toAminoAcidFragment()).collect(
                Collectors.toList());
    }

    private List<AminoAcidFragment> process(final List<Integer> boundaries, final List<NucleotideFragment> fragments)
    {
        /*
                if (fragments.isEmpty()) {
            return listOf()
        }

        val qualEnriched = nucleotideQualEnrichment.enrich(fragments)
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minEvidence, boundaries.filter { it <= MAX_AMINO_ACID_BOUNDARY }.toSet())
        val spliceEnriched = spliceEnricher.enrich(qualEnriched)
        val result = aminoAcidEnricher.enrich(spliceEnriched)

        return result

         */

        if(fragments.isEmpty())
            return Lists.newArrayList();

        List<NucleotideFragment> qualEnriched = mNucleotideQualEnrichment.enrich(fragments);

        NucleotideSpliceEnrichment spliceEnricher = new NucleotideSpliceEnrichment(
                mMinBaseQuality, mMinEvidence, boundaries.stream().filter(x -> x <= MAX_AMINO_ACID_BOUNDARY).collect(Collectors.toSet()));

        List<NucleotideFragment> spliceEnriched = spliceEnricher.enrich(qualEnriched);
        return mAminoAcidEnricher.enrich(spliceEnriched);
    }
}
