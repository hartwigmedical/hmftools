package com.hartwig.hmftools.orange.algo.purple;

import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.jetbrains.annotations.NotNull;

class GermlineGainDeletionUtil
{
    public static double getSomaticMaxCopyNumber(@NotNull List<GermlineAmpDel> gainDelsForGene,
            @NotNull GeneCopyNumber somaticGeneCopyNumber, @NotNull TranscriptData transcript)
    {
        double maximumTumorCopyNumberFromDeletions = getMaximumTumorCopyNumberFromAmpDels(gainDelsForGene);
        double maxCopyNumber = deletionsCoverTranscript(gainDelsForGene, transcript)
                ? maximumTumorCopyNumberFromDeletions
                : Math.max(maximumTumorCopyNumberFromDeletions, somaticGeneCopyNumber.maxCopyNumber());
        return Math.max(0, maxCopyNumber);
    }

    public static double getSomaticMinCopyNumber(@NotNull List<GermlineAmpDel> gainDelsForGene)
    {
        double minCopyNumberFromDeletions = getMinimumTumorCopyNumberFromAmpDels(gainDelsForGene);
        return Math.max(0, minCopyNumberFromDeletions);
    }

    public static double getGermlineMinCopyNumber(@NotNull List<GermlineAmpDel> gainDelsForGene)
    {
        return getCopyNumberFromAmpDels(gainDelsForGene,
                germlineAmpDel -> germlineAmpDel.GermlineCopyNumber,
                Double::min);
    }

    public static boolean deletionsCoverTranscript(@NotNull List<GermlineAmpDel> gainDelsForGene, @NotNull TranscriptData transcript)
    {
        for(ExonData exon : transcript.exons())
        {
            if(!deletionsCoverExon(gainDelsForGene, exon))
            {
                return false;
            }
        }
        return true;
    }

    public static GermlineStatus getGermlineGainDelStatus(@NotNull List<GermlineAmpDel> gainDelsForGene)
    {
        if(gainDelsForGene.isEmpty())
        {
            throw new IllegalArgumentException("Cannot determine germline status for empty list of germline gain dels");
        }

        return gainDelsForGene.stream().map(d -> d.NormalStatus).min(Comparator.naturalOrder()).orElseThrow();
    }

    public static GermlineStatus getSomaticGainDelStatus(@NotNull List<GermlineAmpDel> gainDelsForGene)
    {
        if(gainDelsForGene.isEmpty())
        {
            throw new IllegalArgumentException("Cannot determine germline status for empty list of germline gain dels");
        }

        return gainDelsForGene.stream().map(d -> d.TumorStatus).min(Comparator.naturalOrder()).orElseThrow();
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
    public static TranscriptData findTranscript(@NotNull String geneNameToFind, @NotNull String transcriptNameToFind, @NotNull EnsemblDataCache ensemblDataCache)
    {
        GeneData gene = ensemblDataCache.getGeneDataByName(geneNameToFind);
        if(gene == null)
        {
            throw new IllegalStateException("Could not find gene in ensembl data cache with name: " + geneNameToFind);
        }

        TranscriptData transcript = ensemblDataCache.getTranscriptData(gene.GeneId, transcriptNameToFind);
        if(transcript == null)
        {
            throw new IllegalStateException(
                    String.format("Could not find transcript in ensembl data cache for gene id: %s transcript name: %s",
                            gene.GeneId, transcriptNameToFind));
        }

        return transcript;
    }

    @NotNull
    public static String getChromosome(@NotNull List<GermlineAmpDel> deletions)
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
    public static String getChromosomeBand(@NotNull List<GermlineAmpDel> deletions)
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
        return allSomaticGeneCopyNumbers.stream()
                .filter(g -> g.geneName().equals(geneNameToFind))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("Could not find gene copy number for gene with name: " + geneNameToFind));
    }

    @NotNull
    public static GeneCopyNumber findGeneCopyNumberForGeneTranscript(
            @NotNull String geneNameToFind,
            @NotNull String transcriptToFind,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        return allSomaticGeneCopyNumbers.stream()
                .filter(g -> g.geneName().equals(geneNameToFind) && g.TransName.equals(transcriptToFind))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException(
                        String.format("Could not find gene copy number for gene: %s transcript: %s", geneNameToFind, transcriptToFind)));
    }

    private static double getMinimumTumorCopyNumberFromAmpDels(@NotNull List<GermlineAmpDel> gainDelsForGene)
    {
        return getCopyNumberFromAmpDels(gainDelsForGene,
                germlineAmpDel -> germlineAmpDel.TumorCopyNumber,
                Double::min);
    }

    private static double getMaximumTumorCopyNumberFromAmpDels(@NotNull List<GermlineAmpDel> gainDelsForGene)
    {
        return getCopyNumberFromAmpDels(gainDelsForGene,
                germlineAmpDel -> germlineAmpDel.TumorCopyNumber,
                Double::max);
    }

    private static double getCopyNumberFromAmpDels(
            @NotNull List<GermlineAmpDel> ampDelsForGene,
            Function<GermlineAmpDel, Double> copyNumberGetter,
            BinaryOperator<Double> minMaxOperator)
    {
        if(ampDelsForGene.isEmpty())
        {
            throw new IllegalArgumentException("Cannot determine copy number for empty list of germline amp dels");
        }

        return Math.max(0, ampDelsForGene.stream().map(copyNumberGetter).reduce(minMaxOperator).orElseThrow());
    }

    private static boolean deletionsCoverExon(@NotNull List<GermlineAmpDel> gainDelsForGene, @NotNull ExonData exon)
    {
        List<GermlineAmpDel> relevantSortedDeletions = gainDelsForGene.stream()
                .filter(d -> overlaps(d, exon))
                .sorted(Comparator.comparingInt(o -> o.RegionStart))
                .toList();

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
        for(GermlineAmpDel deletion : relevantSortedDeletions)
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

    private static boolean overlaps(@NotNull GermlineAmpDel deletion, @NotNull ExonData exon)
    {
        return Math.max(deletion.RegionStart, exon.Start) <= Math.min(deletion.RegionEnd, exon.End);
    }
}
