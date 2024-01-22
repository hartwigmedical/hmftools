package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;

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
    public Map<PurpleGainLoss, GermlineDeletion> mapDeletions(@NotNull List<GermlineDeletion> germlineDeletions,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        Map<PurpleGainLoss, GermlineDeletion> deletionMap = Maps.newHashMap();
        for(GermlineDeletion germlineDeletion : germlineDeletions)
        {
            if(germlineDeletion.TumorStatus == GermlineStatus.HOM_DELETION)
            {
                PurpleGainLoss gainLoss = toGainLoss(germlineDeletion, allSomaticGeneCopyNumbers);
                if(!deletionMap.containsKey(gainLoss))
                {
                    deletionMap.put(gainLoss, germlineDeletion);
                }
            }
        }
        return deletionMap;
    }

    @NotNull
    private PurpleGainLoss toGainLoss(@NotNull GermlineDeletion deletion, @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        TranscriptData canonicalTranscript = GermlineCopyNumberUtil.findCanonicalTranscript(deletion.GeneName, ensemblDataCache);
        CopyNumberInterpretation interpretation = GermlineCopyNumberUtil.deletionCoversTranscript(deletion, canonicalTranscript)
                ? CopyNumberInterpretation.FULL_LOSS
                : CopyNumberInterpretation.PARTIAL_LOSS;
        double minCopies = GermlineCopyNumberUtil.getSomaticMinCopyNumber(deletion);
        double maxCopies = GermlineCopyNumberUtil.getSomaticMaxCopyNumber(deletion, allSomaticGeneCopyNumbers, ensemblDataCache);

        return ImmutablePurpleGainLoss.builder()
                .interpretation(interpretation)
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
