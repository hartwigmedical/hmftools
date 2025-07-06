package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public final class Candidates
{
    private final LilacConfig mConfig;
    private final List<HlaSequenceLoci> mNucleotideSequences;
    private final List<HlaSequenceLoci> mAminoAcidSequences;

    public Candidates(final LilacConfig config, final List<HlaSequenceLoci> nucleotideSequences,
            final List<HlaSequenceLoci> aminoAcidSequences)
    {
        mConfig = config;
        mNucleotideSequences = nucleotideSequences;
        mAminoAcidSequences = aminoAcidSequences;
    }

    public List<HlaAllele> unphasedCandidates(final HlaContext context, final Collection<Fragment> fragments,
            final Collection<HlaAllele> commonAllles)
    {
        Set<Integer> aminoAcidBoundary = context.AminoAcidBoundaries;

        List<Fragment> contextFragments = fragments.stream().filter(x -> x.readGene().equals(context.geneName())).toList();
        LL_LOGGER.debug("gene({}) determining un-phased candidates from frags({})", context.geneName(), contextFragments.size());

        SequenceCount aminoAcidCounts = SequenceCount.aminoAcids(mConfig.MinVafFilterDepth, mConfig.MinEvidenceFactor, contextFragments);

        List<HlaSequenceLoci> geneCandidates = mAminoAcidSequences.stream()
                .filter(x -> x.Allele.Gene.equals(context.Gene)).collect(Collectors.toList());

        LL_LOGGER.debug("gene({}) {} candidates before filtering", context.geneName(), geneCandidates.size());

        // Amino acid filtering
        List<HlaSequenceLoci> aminoAcidCandidates = filterSequencesByMinSupport(
	        geneCandidates, aminoAcidCounts, context.AminoAcidBoundaries);

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
        NucleotideFiltering nucleotideFiltering = new NucleotideFiltering(
                mConfig.MinVafFilterDepth, mConfig.MinEvidenceFactor, aminoAcidBoundary);

        List<HlaSequenceLoci> nucleotideCandidatesAfterAminoAcidFiltering = mNucleotideSequences.stream()
                .filter(x -> aminoAcidSpecificAllelesCandidates.contains(x.Allele.asFourDigit()))
                .collect(Collectors.toList());

        List<HlaSequenceLoci> nucleotideSpecificSequences = nucleotideFiltering.filterCandidatesOnAminoAcidBoundaries(
                nucleotideCandidatesAfterAminoAcidFiltering, contextFragments);

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

    private List<HlaSequenceLoci> filterSequencesByMinSupport(
            final Collection<HlaSequenceLoci> candidates, final SequenceCount aminoAcidCount, final Set<Integer> aminoAcidBoundaries)
    {
        // eliminate sequences without min support for their amino acid at each loco, ignoring exon boundaries
        List<HlaSequenceLoci> candidateSequences = Lists.newArrayList();
        candidateSequences.addAll(candidates);

        for(int locus : aminoAcidCount.seqCountsByLoci().keySet())
        {
            if(aminoAcidBoundaries.contains(locus))
            {
                continue;
            }

            List<String> expectedSequences = aminoAcidCount.getMinEvidenceSequences(locus, mConfig.MinEvidenceFactor);

            if(expectedSequences.isEmpty())
            {
                continue;
            }

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
            final HlaContext context, final Collection<HlaAllele> unphasedCandidateAlleles, final Iterable<PhasedEvidence> phasedEvidence)
    {
        LL_LOGGER.debug("gene({}) determining phased candidate set", context.geneName());

        List<HlaSequenceLoci> unphasedCandidates = mAminoAcidSequences.stream()
                .filter(x -> unphasedCandidateAlleles.contains(x.Allele.asFourDigit())).collect(Collectors.toList());

        List<HlaSequenceLoci> phasedCandidates = filterCandidates(unphasedCandidates, phasedEvidence);
        List<HlaAllele> phasedAlleles = phasedCandidates.stream().map(x -> x.Allele).collect(Collectors.toList());

        LL_LOGGER.info("gene({}) has {} candidates after phasing: {}",
                context.geneName(), phasedCandidates.size(), HlaAllele.toString(phasedAlleles, 100));

        return phasedAlleles;
    }

    private static List<HlaSequenceLoci> filterCandidates(final Collection<HlaSequenceLoci> initialCandidates,
            final Iterable<PhasedEvidence> evidence)
    {
        List<HlaSequenceLoci> candidates = Lists.newArrayList();
        candidates.addAll(initialCandidates);

        for(PhasedEvidence newEvidence : evidence)
        {
            candidates = candidates.stream()
                    .filter(x -> x.consistentWithAny(
                            new ArrayList<>(newEvidence.getEvidence().keySet()), newEvidence.getAminoAcidLoci()))
                    .collect(Collectors.toList());
        }

        return candidates;
    }
}
