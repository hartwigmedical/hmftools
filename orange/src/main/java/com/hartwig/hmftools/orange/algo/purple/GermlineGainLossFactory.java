package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;

import org.jetbrains.annotations.NotNull;

public class GermlineGainLossFactory
{
    @NotNull
    private final EnsemblDataCache ensemblDataCache;

    public GermlineGainLossFactory(@NotNull final EnsemblDataCache ensemblDataCache)
    {
        this.ensemblDataCache = ensemblDataCache;
    }

    @NotNull
    public Map<PurpleGainLoss, Boolean> getReportabilityMap(@NotNull List<GermlineDeletion> germlineDeletions,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        List<GermlineDeletion> germlineDeletionsHomozygousInTumor =
                germlineDeletions.stream().filter(d -> d.TumorStatus == GermlineStatus.HOM_DELETION).collect(Collectors.toList());
        Set<String> relevantGeneNames = germlineDeletionsHomozygousInTumor.stream().map(d -> d.GeneName).collect(Collectors.toSet());

        Map<PurpleGainLoss, Boolean> lossToReportability = Maps.newHashMap();
        for(String geneName : relevantGeneNames)
        {
            List<GermlineDeletion> deletionsForGene =
                    germlineDeletionsHomozygousInTumor.stream().filter(d -> d.GeneName.equals(geneName)).collect(Collectors.toList());
            GeneCopyNumber somaticGeneCopyNumber = GermlineDeletionUtil.findGeneCopyNumberForGene(geneName, allSomaticGeneCopyNumbers);

            PurpleGainLoss loss = toGainLoss(geneName, deletionsForGene, somaticGeneCopyNumber);
            boolean reported = deletionsForGene.stream().anyMatch(d -> d.Reported);
            lossToReportability.put(loss, reported);
        }
        return lossToReportability;
    }

    @NotNull
    private PurpleGainLoss toGainLoss(@NotNull String geneName, @NotNull List<GermlineDeletion> deletionsForGene,
            @NotNull GeneCopyNumber somaticGeneCopyNumber)
    {
        TranscriptData canonicalTranscript = GermlineDeletionUtil.findCanonicalTranscript(geneName, ensemblDataCache);
        CopyNumberInterpretation interpretation = GermlineDeletionUtil.deletionsCoverTranscript(deletionsForGene, canonicalTranscript)
                ? CopyNumberInterpretation.FULL_LOSS
                : CopyNumberInterpretation.PARTIAL_LOSS;
        double minCopies = GermlineDeletionUtil.getSomaticMinCopyNumber(deletionsForGene);
        double maxCopies = GermlineDeletionUtil.getSomaticMaxCopyNumber(deletionsForGene, somaticGeneCopyNumber, canonicalTranscript);
        String chromosome = GermlineDeletionUtil.getChromosome(deletionsForGene);
        String chromosomeBand = GermlineDeletionUtil.getChromosomeBand(deletionsForGene);

        return ImmutablePurpleGainLoss.builder()
                .interpretation(interpretation)
                .chromosome(chromosome)
                .chromosomeBand(chromosomeBand)
                .gene(geneName)
                .transcript(canonicalTranscript.TransName)
                .isCanonical(true)
                .minCopies(minCopies)
                .maxCopies(maxCopies)
                .build();
    }
}
