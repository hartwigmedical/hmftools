package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_DEPTH_FILTER;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline.RAW_REF_AMINO_ACID_COUNTS;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public final class Candidates
{
    private final List<HlaSequenceLoci> mNucleotideSequences;
    private final List<HlaSequenceLoci> mAminoAcidSequences;

    public Candidates(final List<HlaSequenceLoci> nucleotideSequences, final List<HlaSequenceLoci> aminoAcidSequences)
    {
        mNucleotideSequences = nucleotideSequences;
        mAminoAcidSequences = aminoAcidSequences;
    }

    public List<HlaAllele> unphasedCandidates(final HlaContext context, final List<Fragment> fragments, final Collection<HlaAllele> commonAllles)
    {
        List<Integer> aminoAcidBoundary = context.AminoAcidBoundaries;

        LL_LOGGER.debug("gene({}) determining un-phased candidates from frags({})", context.geneName(), fragments.size());

        SequenceCount aminoAcidCounts = SequenceCount.buildFromAminoAcids(MIN_EVIDENCE_FACTOR, fragments);

        List<HlaSequenceLoci> geneCandidates = mAminoAcidSequences.stream()
                .filter(x -> x.Allele.Gene.equals(context.Gene)).collect(Collectors.toList());

        LL_LOGGER.debug("gene({}) {} candidates before filtering", context.geneName(), geneCandidates.size());

        // Amino acid filtering
        List<HlaSequenceLoci> aminoAcidCandidates = filterSequencesByMinSupport(geneCandidates, aminoAcidCounts,
                Sets.newTreeSet(context.AminoAcidBoundaries), RAW_REF_AMINO_ACID_COUNTS.get(context.geneName()));

        List<HlaAllele> aminoAcidCandidateAlleles = aminoAcidCandidates.stream().map(x -> x.Allele).collect(Collectors.toList());

        List<HlaAllele> aminoAcidSpecificAllelesCandidates = aminoAcidCandidateAlleles.stream()
                .map(HlaAllele::asFourDigit).toList();

        if(aminoAcidSpecificAllelesCandidates.isEmpty())
        {
            LL_LOGGER.warn("gene({}) no candidates after amino acid filtering - reverting to common allele gene candidates",
                    context.geneName());

            return commonAllles.stream().filter(x -> x.Gene.equals(context.Gene)).collect(Collectors.toList());
        }

        LL_LOGGER.info("gene({}) {} candidates after amino acid filtering", context.geneName(), aminoAcidCandidates.size());

        // Nucleotide filtering
        NucleotideFiltering nucleotideFiltering = new NucleotideFiltering(aminoAcidBoundary);

        List<HlaSequenceLoci> nucleotideCandidatesAfterAminoAcidFiltering = mNucleotideSequences.stream()
                .filter(x -> aminoAcidSpecificAllelesCandidates.contains(x.Allele.asFourDigit()))
                .collect(Collectors.toList());

        List<HlaSequenceLoci> nucleotideSpecificSequences = nucleotideFiltering.filterCandidatesOnAminoAcidBoundaries(
                nucleotideCandidatesAfterAminoAcidFiltering, fragments);

        List<HlaAllele> nucleotideSpecificAllelesCandidates = HlaAllele.dedup
                (nucleotideSpecificSequences.stream().map(x -> x.Allele.asFourDigit()).collect(Collectors.toList()));

        if(nucleotideSpecificAllelesCandidates.isEmpty())
        {
            LL_LOGGER.warn("gene({}) 0 candidates after exon boundary filtering - reverting to amino acid candidates", context.geneName());
            return aminoAcidCandidateAlleles;
        }

        LL_LOGGER.info("gene({}) {} candidates after exon boundary filtering",
                context.geneName(), nucleotideSpecificAllelesCandidates.size());

        return nucleotideSpecificAllelesCandidates;
    }

    @VisibleForTesting
    public static List<HlaSequenceLoci> filterSequencesByMinSupport(final Collection<HlaSequenceLoci> candidates,
            final SequenceCount aminoAcidCount, final Set<Integer> aminoAcidBoundaries, final SequenceCount rawAminoAcidCounts)
    {
        // eliminate sequences without min support for their amino acid at each loco, ignoring exon boundaries
        List<HlaSequenceLoci> candidateSequences = Lists.newArrayList();
        candidateSequences.addAll(candidates);

        for(final int locus : aminoAcidCount.seqCountsByLoci().keySet())
        {
            if(aminoAcidBoundaries.contains(locus))
                continue;

            if(rawAminoAcidCounts.get(locus).size() < MIN_DEPTH_FILTER)
                continue;

            Set<String> expectedSequences = Sets.newHashSet(aminoAcidCount.getMinEvidenceSequences(locus));
            if(expectedSequences.isEmpty())
                continue;

            int index = 0;
            while(index < candidateSequences.size())
            {
                HlaSequenceLoci sequence = candidateSequences.get(index);

                if(locus < sequence.length())
                {
                    String locusSeq = sequence.sequence(locus);

                    if(locusSeq.equals(WILD_STR) || expectedSequences.contains(locusSeq))
                    {
                        ++index;
                        continue;
                    }
                }

                candidateSequences.remove(index);
            }
        }

        return candidateSequences;
    }

    public List<HlaAllele> phasedCandidates(
            final HlaContext context, final Set<HlaAllele> unphasedCandidateAlleles, final Iterable<PhasedEvidence> phasedEvidence)
    {
        LL_LOGGER.debug("gene({}) determining phased candidate set", context.geneName());

        List<HlaSequenceLoci> unphasedCandidates = mAminoAcidSequences.stream()
                .filter(x -> unphasedCandidateAlleles.contains(x.Allele.asFourDigit())).collect(Collectors.toList());

        List<HlaSequenceLoci> phasedCandidates = filterCandidates(
                unphasedCandidates, phasedEvidence, RAW_REF_AMINO_ACID_COUNTS.get(context.geneName()));
        List<HlaAllele> phasedAlleles = phasedCandidates.stream().map(x -> x.Allele).collect(Collectors.toList());

        LL_LOGGER.info("gene({}) has {} candidates after phasing: {}",
                context.geneName(), phasedCandidates.size(), HlaAllele.toString(phasedAlleles, 100));

        return phasedAlleles;
    }

    @VisibleForTesting
    public static List<HlaSequenceLoci> filterCandidates(
            final Collection<HlaSequenceLoci> initialCandidates, final Iterable<PhasedEvidence> evidence, final SequenceCount rawAminoAcidCounts)
    {
        List<HlaSequenceLoci> candidates = Lists.newArrayList();
        candidates.addAll(initialCandidates);

        for(PhasedEvidence newEvidence : evidence)
        {
            List<Integer> targetLoci = newEvidence.getAminoAcidLoci();
            List<StringBuilder> targetSequenceBuilders = newEvidence.getEvidence().keySet().stream()
                    .map(StringBuilder::new)
                    .toList();

            List<Integer> lowDepthIndices = Lists.newArrayList();
            for(int j = 0; j < targetLoci.size(); j++)
            {
                int locus = targetLoci.get(j);
                if(rawAminoAcidCounts.get(locus).size() < MIN_DEPTH_FILTER)
                    lowDepthIndices.add(j);
            }

            if(lowDepthIndices.size() >= targetLoci.size() - 1)
                continue;

            if(!lowDepthIndices.isEmpty())
            {
                Collections.sort(lowDepthIndices, Comparator.comparingInt((Integer x) -> x).reversed());
                for(int index : lowDepthIndices)
                {
                    targetLoci.remove(index);
                    for(StringBuilder seqBuilder : targetSequenceBuilders)
                        seqBuilder.deleteCharAt(index);
                }
            }

            List<String> targetSequences = targetSequenceBuilders.stream().map(StringBuilder::toString).collect(Collectors.toList());
            candidates = candidates.stream().filter(x -> x.consistentWithAny(targetSequences, targetLoci)).collect(Collectors.toList());
        }

        return candidates;
    }
}
