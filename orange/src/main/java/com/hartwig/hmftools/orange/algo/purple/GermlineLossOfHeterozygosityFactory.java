package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;

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
    public Map<PurpleLossOfHeterozygosity, GermlineDeletion> mapDeletions(@NotNull List<GermlineDeletion> germlineDeletions,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        Map<PurpleLossOfHeterozygosity, GermlineDeletion> deletionMap = Maps.newHashMap();
        for(GermlineDeletion germlineDeletion : germlineDeletions)
        {
            if(germlineDeletion.TumorStatus == GermlineStatus.HET_DELETION)
            {
                PurpleLossOfHeterozygosity lossOfHeterozygosity = toLossOfHeterozygosity(germlineDeletion, allSomaticGeneCopyNumbers);
                if(!deletionMap.containsKey(lossOfHeterozygosity))
                {
                    deletionMap.put(lossOfHeterozygosity, germlineDeletion);
                }
            }
        }
        return deletionMap;
    }

    @NotNull
    private PurpleLossOfHeterozygosity toLossOfHeterozygosity(@NotNull GermlineDeletion deletion,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        TranscriptData canonicalTranscript = GermlineCopyNumberUtil.findCanonicalTranscript(deletion.GeneName, ensemblDataCache);
        GeneProportion geneProportion = GermlineCopyNumberUtil.deletionCoversTranscript(deletion, canonicalTranscript)
                ? GeneProportion.FULL_GENE
                : GeneProportion.PARTIAL_GENE;
        double minCopies = GermlineCopyNumberUtil.getSomaticMinCopyNumber(deletion);
        double maxCopies = GermlineCopyNumberUtil.getSomaticMaxCopyNumber(deletion, allSomaticGeneCopyNumbers, ensemblDataCache);

        return ImmutablePurpleLossOfHeterozygosity.builder()
                .geneProportion(geneProportion)
                .chromosome(deletion.Chromosome)
                .chromosomeBand(deletion.ChromosomeBand)
                .gene(deletion.GeneName)
                .transcript(canonicalTranscript.TransName)
                .isCanonical(true)
                .minCopies(minCopies)
                .maxCopies(maxCopies)
                .build();
    }
}
