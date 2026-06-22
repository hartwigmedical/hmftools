package com.hartwig.hmftools.tars.liftback.rescue;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.tars.common.ContigEntry;

// Derives annotated splice junctions from the ensembl data cache into a Set<ChrBaseRegion> for O(1) lookup.
// Pairs adjacent exons of each transcript to derive intron coords; duplicates across transcripts collapse in the Set.
public final class AnnotatedJunctionLoader
{
    private AnnotatedJunctionLoader() { }

    public static Set<ChrBaseRegion> load(final String ensemblDir, final RefGenomeVersion refGenomeVersion)
    {
        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, false);
        ensemblDataCache.load(false);
        return deriveIntrons(ensemblDataCache, refGenomeVersion);
    }

    // derive the same intron set from the contig sidecar: each entry's exonSpans are a multi-exon transcript's
    // genomic exons. Single-exon transcripts contribute no junctions, so the sidecar misses nothing here.
    public static Set<ChrBaseRegion> fromContigEntries(final List<ContigEntry> entries)
    {
        Set<ChrBaseRegion> introns = new HashSet<>();
        for(final ContigEntry entry : entries)
        {
            addSpanIntrons(introns, entry.chromosome(), entry.exonSpans());
        }
        return introns;
    }

    private static void addSpanIntrons(final Set<ChrBaseRegion> introns, final String chromosome, final List<BaseRegion> exonSpans)
    {
        if(exonSpans.size() < 2)
        {
            return;
        }

        List<BaseRegion> exons = new ArrayList<>(exonSpans);
        exons.sort(Comparator.comparingInt(BaseRegion::start));

        for(int i = 0; i < exons.size() - 1; ++i)
        {
            int intronStart = exons.get(i).end() + 1;
            int intronEnd = exons.get(i + 1).start() - 1;
            if(intronEnd < intronStart)
                continue;     // abutting exons -> no intron

            introns.add(new ChrBaseRegion(chromosome, intronStart, intronEnd));
        }
    }

    public static Set<ChrBaseRegion> deriveIntrons(final EnsemblDataCache ensemblDataCache, final RefGenomeVersion refGenomeVersion)
    {
        Set<ChrBaseRegion> introns = new HashSet<>();
        for(final List<GeneData> genes : ensemblDataCache.getChrGeneDataMap().values())
        {
            for(final GeneData gene : genes)
            {
                List<TranscriptData> transcripts = ensemblDataCache.getTranscripts(gene.GeneId);
                if(transcripts == null)
                    continue;

                String chromosome = refGenomeVersion.versionedChromosome(gene.Chromosome);
                for(final TranscriptData transcript : transcripts)
                {
                    addTranscriptIntrons(introns, chromosome, transcript.exons());
                }
            }
        }
        return introns;
    }

    private static void addTranscriptIntrons(final Set<ChrBaseRegion> introns, final String chromosome,
            final List<ExonData> transcriptExons)
    {
        if(transcriptExons.size() < 2)
        {
            return;
        }

        List<ExonData> exons = new ArrayList<>(transcriptExons);
        exons.sort(Comparator.comparingInt(exon -> exon.Rank));

        for(int i = 0; i < exons.size() - 1; ++i)
        {
            ExonData first = exons.get(i);
            ExonData second = exons.get(i + 1);
            // Rank follows transcription direction, not genomic coords -- derive intron bounds from coords directly.
            int intronStart = Math.min(first.End, second.End) + 1;
            int intronEnd = Math.max(first.Start, second.Start) - 1;
            if(intronEnd < intronStart)
                continue;     // abutting or overlapping exons -> no intron

            introns.add(new ChrBaseRegion(chromosome, intronStart, intronEnd));
        }
    }
}
