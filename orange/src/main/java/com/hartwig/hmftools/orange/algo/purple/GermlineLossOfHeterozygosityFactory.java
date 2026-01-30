package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.GeneProportion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;

import org.jetbrains.annotations.NotNull;

public class GermlineLossOfHeterozygosityFactory
{
    private final EnsemblDataCache ensemblDataCache;

    public GermlineLossOfHeterozygosityFactory(@NotNull final EnsemblDataCache ensemblDataCache)
    {
        this.ensemblDataCache = ensemblDataCache;
    }

    public Map<PurpleLossOfHeterozygosity, Boolean> getReportabilityMap(@NotNull List<GermlineAmpDel> germlineDeletions,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        List<GermlineAmpDel> germlineDeletionsHeterozygousInTumor =
                germlineDeletions.stream().filter(d -> d.TumorStatus == GermlineStatus.HET_DELETION).collect(Collectors.toList());
        Set<String> relevantGeneNames = germlineDeletionsHeterozygousInTumor.stream().map(d -> d.GeneName).collect(Collectors.toSet());

        Map<PurpleLossOfHeterozygosity, Boolean> lohToReportability = Maps.newHashMap();
        for(String geneName : relevantGeneNames)
        {
            List<GermlineAmpDel> deletionsForGene =
                    germlineDeletionsHeterozygousInTumor.stream().filter(d -> d.GeneName.equals(geneName)).collect(Collectors.toList());
            GeneCopyNumber somaticGeneCopyNumber = GermlineGainDeletionUtil.findGeneCopyNumberForGene(geneName, allSomaticGeneCopyNumbers);

            PurpleLossOfHeterozygosity  lossOfHeterozygosity = toPurpleLossOfHeterozygosity(
                    geneName, deletionsForGene, somaticGeneCopyNumber);

            boolean reported = deletionsForGene.stream().anyMatch(d -> d.Reported == ReportedStatus.REPORTED);
            lohToReportability.put(lossOfHeterozygosity, reported);
        }
        return lohToReportability;
    }

    @NotNull
    private PurpleLossOfHeterozygosity toPurpleLossOfHeterozygosity(@NotNull String geneName,
            @NotNull List<GermlineAmpDel> deletionsForGene, @NotNull GeneCopyNumber somaticGeneCopyNumber)
    {
        TranscriptData canonicalTranscript = GermlineGainDeletionUtil.findCanonicalTranscript(geneName, ensemblDataCache);
        GeneProportion geneProportion = GermlineGainDeletionUtil.deletionsCoverTranscript(deletionsForGene, canonicalTranscript)
                ? GeneProportion.FULL_GENE
                : GeneProportion.PARTIAL_GENE;
        double minCopies = GermlineGainDeletionUtil.getSomaticMinCopyNumber(deletionsForGene);
        double maxCopies = GermlineGainDeletionUtil.getSomaticMaxCopyNumber(deletionsForGene, somaticGeneCopyNumber, canonicalTranscript);
        String chromosome = GermlineGainDeletionUtil.getChromosome(deletionsForGene);
        String chromosomeBand = GermlineGainDeletionUtil.getChromosomeBand(deletionsForGene);

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
