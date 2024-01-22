package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;

import org.jetbrains.annotations.NotNull;

public class GermlineCopyNumberUtil
{

    public static double getSomaticMaxCopyNumber(@NotNull GermlineDeletion deletion,
            @NotNull List<GeneCopyNumber> allSomaticGeneCopyNumbers, EnsemblDataCache ensemblDataCache)
    {
        TranscriptData canonicalTranscript = findCanonicalTranscript(deletion.GeneName, ensemblDataCache);
        return deletionCoversTranscript(deletion, canonicalTranscript)
                ? getSomaticMinCopyNumber(deletion)
                : maxCopyNumberFromGeneCopyNumber(deletion.GeneName, allSomaticGeneCopyNumbers);
    }

    public static double getSomaticMinCopyNumber(@NotNull GermlineDeletion deletion)
    {
        return Math.max(0, deletion.TumorCopyNumber);
    }

    public static boolean deletionCoversTranscript(@NotNull GermlineDeletion deletion, @NotNull TranscriptData canonicalTranscript)
    {
        return deletion.RegionStart < canonicalTranscript.TransStart && deletion.RegionEnd > canonicalTranscript.TransEnd;
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
