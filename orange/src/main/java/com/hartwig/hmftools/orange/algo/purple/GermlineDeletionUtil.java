package com.hartwig.hmftools.orange.algo.purple;

import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;

import org.jetbrains.annotations.NotNull;

class GermlineDeletionUtil
{
    public static double getSomaticMaxCopyNumber(@NotNull List<GermlineDeletion> deletionsForGene,
            @NotNull GeneCopyNumber somaticGeneCopyNumber, @NotNull TranscriptData transcript)
    {
        double maximumTumorCopyNumberFromDeletions = getMaximumTumorCopyNumberFromDeletions(deletionsForGene);
        double maxCopyNumber = deletionsCoverTranscript(deletionsForGene, transcript)
                ? maximumTumorCopyNumberFromDeletions
                : Math.max(maximumTumorCopyNumberFromDeletions, somaticGeneCopyNumber.maxCopyNumber());
        return Math.max(0, maxCopyNumber);
    }

    public static double getSomaticMinCopyNumber(@NotNull List<GermlineDeletion> deletionsForGene)
    {
        double minCopyNumberFromDeletions = getMinimumTumorCopyNumberFromDeletions(deletionsForGene);
        return Math.max(0, minCopyNumberFromDeletions);
    }

    public static boolean deletionsCoverTranscript(@NotNull List<GermlineDeletion> deletionsForGene, @NotNull TranscriptData transcript)
    {
        for(ExonData exon : transcript.exons())
        {
            if(!deletionsCoverExon(deletionsForGene, exon))
            {
                return false;
            }
        }
        return true;
    }

    @NotNull
    public static TranscriptData findCanonicalTranscript(@NotNull String geneNameToFind, @NotNull EnsemblDataCache ensemblDataCache)
    {
        GeneData gene = ensemblDataCache.getGeneDataByName(geneNameToFind);
        if(gene == null)
        {
            throw new IllegalStateException("Could not find gene in ensembl data cache with name: " + geneNameToFind);
        }

        TranscriptData transcript = ensemblDataCache.getCanonicalTranscriptData(gene.GeneId);
        if(transcript == null)
        {
            throw new IllegalStateException("Could not find canonical transcript in ensembl data cache for gene with id: " + gene.GeneId);
        }

        return transcript;
    }

    @NotNull
    public static String getChromosome(@NotNull List<GermlineDeletion> deletions)
    {
        if(deletions.isEmpty())
        {
            throw new IllegalArgumentException("No germline deletions provided to get chromosome from.");
        }

        Set<String> chromosomes = deletions.stream().map(d -> d.Chromosome).collect(Collectors.toSet());
        if(chromosomes.size() == 1)
        {
            return chromosomes.iterator().next();
        }
        else
        {
            throw new IllegalStateException("Germline deletions disagree on chromosome.");
        }
    }

    @NotNull
    public static String getChromosomeBand(@NotNull List<GermlineDeletion> deletions)
    {
        if(deletions.isEmpty())
        {
            throw new IllegalArgumentException("No germline deletions provided to get chromosome band from.");
        }

        Set<String> chromosomeBands = deletions.stream().map(d -> d.ChromosomeBand).collect(Collectors.toSet());
        if(chromosomeBands.size() == 1)
        {
            return chromosomeBands.iterator().next();
        }
        else
        {
            throw new IllegalStateException("Germline deletions disagree on chromosome band.");
        }
    }

    @NotNull
    public static GeneCopyNumber findGeneCopyNumberForGene(@NotNull String geneNameToFind,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        for(GeneCopyNumber geneCopyNumber : allSomaticGeneCopyNumbers)
        {
            if(geneCopyNumber.geneName().equals(geneNameToFind))
            {
                return geneCopyNumber;
            }
        }

        throw new IllegalStateException("Could not find gene copy number for gene with name: " + geneNameToFind);
    }

    private static double getMinimumTumorCopyNumberFromDeletions(@NotNull List<GermlineDeletion> deletionsForGene)
    {
        if(deletionsForGene.isEmpty())
        {
            throw new IllegalArgumentException("Cannot determine minimum tumor copy number for empty list of germline deletions");
        }

        return deletionsForGene.stream().mapToDouble(d -> d.TumorCopyNumber).min().getAsDouble();
    }

    private static double getMaximumTumorCopyNumberFromDeletions(@NotNull List<GermlineDeletion> deletionsForGene)
    {
        if(deletionsForGene.isEmpty())
        {
            throw new IllegalArgumentException("Cannot determine maximum tumor copy number for empty list of germline deletions");
        }

        return deletionsForGene.stream().mapToDouble(d -> d.TumorCopyNumber).max().getAsDouble();
    }

    private static boolean deletionsCoverExon(@NotNull List<GermlineDeletion> deletionsForGene, @NotNull ExonData exon)
    {
        List<GermlineDeletion> relevantSortedDeletions = deletionsForGene.stream()
                .filter(d -> overlaps(d, exon))
                .sorted(Comparator.comparingInt(o -> o.RegionStart))
                .collect(Collectors.toList());

        if(relevantSortedDeletions.isEmpty())
        {
            return false;
        }

        boolean startPositionCovered = relevantSortedDeletions.get(0).RegionStart <= exon.Start;
        if(!startPositionCovered)
        {
            return false;
        }

        int latestPositionCoveredContinuously = relevantSortedDeletions.get(0).RegionEnd;
        for(GermlineDeletion deletion : relevantSortedDeletions)
        {
            if(deletion.RegionStart > latestPositionCoveredContinuously + 1)
            {
                return false;
            }
            latestPositionCoveredContinuously = Math.max(latestPositionCoveredContinuously, deletion.RegionEnd);
            if(latestPositionCoveredContinuously >= exon.End)
            {
                return true;
            }
        }

        return false;
    }

    private static boolean overlaps(@NotNull GermlineDeletion deletion, @NotNull ExonData exon)
    {
        return Math.max(deletion.RegionStart, exon.Start) <= Math.min(deletion.RegionEnd, exon.End);
    }
}
