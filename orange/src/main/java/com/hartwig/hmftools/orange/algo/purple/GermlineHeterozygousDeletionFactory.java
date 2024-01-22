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
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleHeterozygousDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleHeterozygousDeletion;

import org.jetbrains.annotations.NotNull;

public class GermlineHeterozygousDeletionFactory
{
    @NotNull
    private final EnsemblDataCache ensemblDataCache;

    public GermlineHeterozygousDeletionFactory(@NotNull final EnsemblDataCache ensemblDataCache)
    {
        this.ensemblDataCache = ensemblDataCache;
    }

    @NotNull
    public Map<PurpleHeterozygousDeletion, GermlineDeletion> mapDeletions(@NotNull List<GermlineDeletion> germlineDeletions,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        Map<PurpleHeterozygousDeletion, GermlineDeletion> deletionMap = Maps.newHashMap();
        for(GermlineDeletion germlineDeletion : germlineDeletions)
        {
            if(germlineDeletion.TumorStatus == GermlineStatus.HET_DELETION)
            {
                PurpleHeterozygousDeletion heterozygousDeletion = toHeterozygousDeletion(germlineDeletion, allSomaticGeneCopyNumbers);
                if(!deletionMap.containsKey(heterozygousDeletion))
                {
                    deletionMap.put(heterozygousDeletion, germlineDeletion);
                }
            }
        }
        return deletionMap;
    }

    @NotNull
    private PurpleHeterozygousDeletion toHeterozygousDeletion(@NotNull GermlineDeletion deletion,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        TranscriptData canonicalTranscript = GermlineCopyNumberUtil.findCanonicalTranscript(deletion.GeneName, ensemblDataCache);
        GeneProportion geneProportion = GermlineCopyNumberUtil.deletionCoversTranscript(deletion, canonicalTranscript)
                ? GeneProportion.FULL_GENE
                : GeneProportion.PARTIAL_GENE;
        double minCopies = GermlineCopyNumberUtil.getSomaticMinCopyNumber(deletion);
        double maxCopies = GermlineCopyNumberUtil.getSomaticMaxCopyNumber(deletion, allSomaticGeneCopyNumbers, ensemblDataCache);

        return ImmutablePurpleHeterozygousDeletion.builder()
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
