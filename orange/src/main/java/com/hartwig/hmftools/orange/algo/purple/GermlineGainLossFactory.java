package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
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
        TranscriptData canonicalTranscript = findCanonicalTranscript(deletion.GeneName);

        boolean isFullDeletion = deletion.RegionStart < canonicalTranscript.TransStart && deletion.RegionEnd > canonicalTranscript.TransEnd;
        double minCopies = Math.max(0, deletion.TumorCopyNumber);
        double maxCopies = isFullDeletion ? minCopies : maxCopyNumberFromGeneCopyNumber(deletion.GeneName, allSomaticGeneCopyNumbers);

        return ImmutablePurpleGainLoss.builder()
                .interpretation(isFullDeletion ? CopyNumberInterpretation.FULL_LOSS : CopyNumberInterpretation.PARTIAL_LOSS)
                .chromosome(deletion.Chromosome)
                .chromosomeBand(deletion.ChromosomeBand)
                .gene(deletion.GeneName)
                .transcript(canonicalTranscript.TransName)
                .isCanonical(true)
                .minCopies(minCopies)
                .maxCopies(maxCopies)
                .build();
    }

    @NotNull
    private TranscriptData findCanonicalTranscript(@NotNull String geneNameToFind)
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

    private static double maxCopyNumberFromGeneCopyNumber(@NotNull String geneNameToFind,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers)
    {
        for(GeneCopyNumber geneCopyNumber : allSomaticGeneCopyNumbers)
        {
            if(geneCopyNumber.geneName().equals(geneNameToFind))
            {
                return Math.max(0, geneCopyNumber.maxCopyNumber());
            }
        }

        throw new IllegalStateException("Could not find gene copy number for gene with name: " + geneNameToFind);
    }
}
