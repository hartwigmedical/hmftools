package com.hartwig.hmftools.lilac.candidates;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragment.nucFragments;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragment;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public final class Candidates
{
    private final LilacConfig mConfig;
    private final List<HlaSequenceLoci> mNucleotideSequences;
    private final List<HlaSequenceLoci> mAminoAcidSequences;

    public Candidates(
            final LilacConfig config, final List<HlaSequenceLoci> nucleotideSequences, final List<HlaSequenceLoci> aminoAcidSequences)
    {
        mConfig = config;
        mNucleotideSequences = nucleotideSequences;
        mAminoAcidSequences = aminoAcidSequences;
    }

    public final List<HlaAllele> unphasedCandidates(final HlaContext context, final List<AminoAcidFragment> fragments)
    {
        List<Integer> aminoAcidBoundary = context.AminoAcidBoundaries;

        LL_LOGGER.info("determining un-phased candidate set for gene {}", context.geneName());

        SequenceCount aminoAcidCounts = SequenceCount.aminoAcids(mConfig.MinEvidence, fragments);

        List<HlaSequenceLoci> geneCandidates = mAminoAcidSequences.stream()
                .filter(x -> x.getAllele().Gene.equals(context.Gene)).collect(Collectors.toList());

        LL_LOGGER.info("  {} candidates before filtering", geneCandidates.size());

        // Amino acid filtering
        List<HlaSequenceLoci> aminoAcidCandidates = filterAminoAcidsByBoundries(geneCandidates, aminoAcidCounts, context.AminoAcidBoundaries);

        List<HlaAllele> aminoAcidCandidateAlleles = aminoAcidCandidates.stream().map(x -> x.getAllele()).collect(Collectors.toList());

        List<HlaAllele> aminoAcidSpecificAllelesCandidates = aminoAcidCandidateAlleles.stream()
                .map(x -> x.asFourDigit()).collect(Collectors.toList());

        if (aminoAcidSpecificAllelesCandidates.isEmpty())
        {
            LL_LOGGER.warn("  0 candidates after amino acid filtering - reverting to all gene candidates");
            return geneCandidates.stream().map(x -> x.getAllele()).collect(Collectors.toList());
        }

        LL_LOGGER.info("  {} candidates after amino acid filtering", aminoAcidCandidates.size());

        // Nucleotide filtering
        NucleotideFiltering nucleotideFiltering = new NucleotideFiltering(mConfig.MinEvidence, aminoAcidBoundary);

        List<HlaSequenceLoci> nucleotideCandidatesAfterAminoAcidFiltering = mNucleotideSequences.stream()
                .filter(x -> aminoAcidSpecificAllelesCandidates.contains(x.getAllele().asFourDigit()))
                .collect(Collectors.toList());

        List<HlaSequenceLoci> nucleotideSpecificSequences = nucleotideFiltering.filterCandidatesOnAminoAcidBoundaries(
                nucleotideCandidatesAfterAminoAcidFiltering, nucFragments(fragments));

        List<HlaAllele> nucleotideSpecificAllelesCandidates = HlaAllele.dedup
                (nucleotideSpecificSequences.stream().map(x -> x.getAllele().asFourDigit()).collect(Collectors.toList()));

        if (nucleotideSpecificAllelesCandidates.isEmpty())
        {
            LL_LOGGER.warn("  0 candidates after exon boundary filtering - reverting to amino acid candidates");
            return aminoAcidCandidateAlleles;
        }

        LL_LOGGER.info("  {} candidates after exon boundary filtering", nucleotideSpecificAllelesCandidates.size());
        return nucleotideSpecificAllelesCandidates;
    }

    private List<HlaSequenceLoci> filterAminoAcidsByBoundries(
            final List<HlaSequenceLoci> candidates, final SequenceCount aminoAcidCount, final List<Integer> aminoAcidBoundaries)
    {
        List<HlaSequenceLoci> results = Lists.newArrayList();
        results.addAll(candidates);

        for(int loci = 0; loci < aminoAcidCount.getLength(); ++loci)
        {
            if(aminoAcidBoundaries.contains(loci))
                continue;

            List<String> expectedSequences = aminoAcidCount.getMinCountSequences(loci);

            final int lociConst = loci;

            results = results.stream()
                    .filter(x -> x.consistentWithAny(expectedSequences, Lists.newArrayList(lociConst)))
                    .collect(Collectors.toList());
        }

        return results;
    }

    public List<HlaAllele> phasedCandidates(
            final HlaContext context, final List<HlaAllele> unphasedCandidateAlleles, final List<PhasedEvidence> phasedEvidence)
    {
        LL_LOGGER.info("determining phased candidate set for gene {}", context.geneName());

        List<HlaSequenceLoci> unphasedCandidates = mAminoAcidSequences.stream()
                .filter(x -> unphasedCandidateAlleles.contains(x.getAllele().asFourDigit())).collect(Collectors.toList());

        List<HlaSequenceLoci> phasedCandidates = filterCandidates(unphasedCandidates, phasedEvidence);
        List<HlaAllele> phasedAlleles = phasedCandidates.stream().map(x -> x.getAllele()).collect(Collectors.toList());

        LL_LOGGER.info("  {} candidates after phasing: {}", phasedCandidates.size(), phasedAlleles);

        return phasedAlleles;
    }

    private List<HlaSequenceLoci> filterCandidates(final List<HlaSequenceLoci> initialCandidates, final List<PhasedEvidence> evidence)
    {
        List<HlaSequenceLoci> candidates = Lists.newArrayList();
        candidates.addAll(initialCandidates);

        for(int i = 0; i < evidence.size(); ++i)
        {
            PhasedEvidence newEvidence = evidence.get(i);

            candidates = candidates.stream()
                .filter(x -> x.consistentWithAny(newEvidence.getEvidence().keySet(), newEvidence.getAminoAcidIndices()))
                .collect(Collectors.toList());
        }

        return candidates;
    }

}
