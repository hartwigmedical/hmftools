package com.hartwig.hmftools.tars.liftback.rescue;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

// Derives annotated splice junctions from the ensembl data cache into a Set<ChrBaseRegion> for O(1) lookup.
// Pairs adjacent exons of each transcript to derive intron coords; duplicates across transcripts collapse in the Set.
public final class AnnotatedJunctionLoader
{
    private AnnotatedJunctionLoader() {}

    public static Set<ChrBaseRegion> load(final String ensemblDir, final RefGenomeVersion refGenomeVersion)
    {
        final EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, false);
        ensemblDataCache.load(false);
        return deriveIntrons(ensemblDataCache, refGenomeVersion);
    }

    public static Set<ChrBaseRegion> deriveIntrons(final EnsemblDataCache ensemblDataCache, final RefGenomeVersion refGenomeVersion)
    {
        final Set<ChrBaseRegion> introns = new HashSet<>();
        for(final List<GeneData> genes : ensemblDataCache.getChrGeneDataMap().values())
        {
            for(final GeneData gene : genes)
            {
                final List<TranscriptData> transcripts = ensemblDataCache.getTranscripts(gene.GeneId);
                if(transcripts == null)
                    continue;

                final String chromosome = refGenomeVersion.versionedChromosome(gene.Chromosome);
                for(final TranscriptData transcript : transcripts)
                    addTranscriptIntrons(introns, chromosome, transcript.exons());
            }
        }
        return introns;
    }

    private static void addTranscriptIntrons(final Set<ChrBaseRegion> introns, final String chromosome, final List<ExonData> transcriptExons)
    {
        if(transcriptExons.size() < 2)
            return;

        final List<ExonData> exons = new ArrayList<>(transcriptExons);
        exons.sort(Comparator.comparingInt(exon -> exon.Rank));

        for(int i = 0; i < exons.size() - 1; ++i)
        {
            final ExonData first = exons.get(i);
            final ExonData second = exons.get(i + 1);
            // Rank follows transcription direction, not genomic coords -- derive intron bounds from coords directly.
            final int intronStart = Math.min(first.End, second.End) + 1;
            final int intronEnd = Math.max(first.Start, second.Start) - 1;
            if(intronEnd < intronStart)
                continue;     // abutting or overlapping exons -> no intron

            introns.add(new ChrBaseRegion(chromosome, intronStart, intronEnd));
        }
    }
}
