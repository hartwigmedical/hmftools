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
import com.hartwig.hmftools.datamodel.purple.GeneProportion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;

import org.jetbrains.annotations.NotNull;

public class GermlineLossOfHeterozygosityFactory
{
    @NotNull
    private final EnsemblDataCache ensemblDataCache;

    public GermlineLossOfHeterozygosityFactory(@NotNull final EnsemblDataCache ensemblDataCache)
    {
        this.ensemblDataCache = ensemblDataCache;
    }

    @NotNull
    public Map<PurpleLossOfHeterozygosity, Boolean> getReportabilityMap(@NotNull List<GermlineDeletion> germlineDeletions,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        List<GermlineDeletion> germlineDeletionsHeterozygousInTumor =
                germlineDeletions.stream().filter(d -> d.TumorStatus == GermlineStatus.HET_DELETION).collect(Collectors.toList());
        Set<String> relevantGeneNames = germlineDeletionsHeterozygousInTumor.stream().map(d -> d.GeneName).collect(Collectors.toSet());

        Map<PurpleLossOfHeterozygosity, Boolean> lohToReportability = Maps.newHashMap();
        for(String geneName : relevantGeneNames)
        {
            List<GermlineDeletion> deletionsForGene =
                    germlineDeletionsHeterozygousInTumor.stream().filter(d -> d.GeneName.equals(geneName)).collect(Collectors.toList());
            GeneCopyNumber somaticGeneCopyNumber = GermlineDeletionUtil.findGeneCopyNumberForGene(geneName, allSomaticGeneCopyNumbers);

            PurpleLossOfHeterozygosity lossOfHeterozygosity =
                    toPurpleLossOfHeterozygosity(geneName, deletionsForGene, somaticGeneCopyNumber);
            boolean reported = deletionsForGene.stream().anyMatch(d -> d.Reported);
            lohToReportability.put(lossOfHeterozygosity, reported);
        }
        return lohToReportability;
    }

    @NotNull
    private PurpleLossOfHeterozygosity toPurpleLossOfHeterozygosity(@NotNull String geneName,
            @NotNull List<GermlineDeletion> deletionsForGene, @NotNull GeneCopyNumber somaticGeneCopyNumber)
    {
        TranscriptData canonicalTranscript = GermlineDeletionUtil.findCanonicalTranscript(geneName, ensemblDataCache);
        GeneProportion geneProportion = GermlineDeletionUtil.deletionsCoverTranscript(deletionsForGene, canonicalTranscript)
                ? GeneProportion.FULL_GENE
                : GeneProportion.PARTIAL_GENE;
        double minCopies = GermlineDeletionUtil.getSomaticMinCopyNumber(deletionsForGene);
        double maxCopies = GermlineDeletionUtil.getSomaticMaxCopyNumber(deletionsForGene, somaticGeneCopyNumber, canonicalTranscript);
        String chromosome = GermlineDeletionUtil.getChromosome(deletionsForGene);
        String chromosomeBand = GermlineDeletionUtil.getChromosomeBand(deletionsForGene);

        return ImmutablePurpleLossOfHeterozygosity.builder()
                .geneProportion(geneProportion)
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
