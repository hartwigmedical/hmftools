package com.hartwig.hmftools.redux.jitter;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.collect.ImmutableIntervalTree;
import com.hartwig.hmftools.common.region.BaseRegion;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMRecord;

public class SampleReadProcessor
{
    private final List<MicrosatelliteSiteData> mMicrosatelliteSiteAnalysers;
    private final Map<String, ImmutableIntervalTree<MicrosatelliteSiteData>> mMicrosatelliteSiteAnalysersByChromosome;

    public SampleReadProcessor(
            final MsJitterConfig config, final Collection<MicrosatelliteSite> microsatelliteSites,
            final ConsensusMarker consensusMarker)
    {
        mMicrosatelliteSiteAnalysers = microsatelliteSites.stream()
                .map(x -> new MicrosatelliteSiteData(x, consensusMarker, config.WriteSiteFile))
                .collect(Collectors.toList());

        mMicrosatelliteSiteAnalysersByChromosome = Maps.newHashMap();
        Multimap<String, MicrosatelliteSiteData> analysersByChromosome =
                Multimaps.index(mMicrosatelliteSiteAnalysers, analyser -> analyser.refGenomeMicrosatellite().chromosome());

        for(Map.Entry<String, Collection<MicrosatelliteSiteData>> chromosomeAnalysers : analysersByChromosome.asMap().entrySet())
        {
            String chromosome = chromosomeAnalysers.getKey();
            Collection<MicrosatelliteSiteData> analysers = chromosomeAnalysers.getValue();

            Collection<Pair<BaseRegion, MicrosatelliteSiteData>> entries = analysers
                    .stream()
                    .map(analyser -> Pair.of(analyser.refGenomeMicrosatellite().Region.baseRegion(), analyser))
                    .collect(Collectors.toList());

            mMicrosatelliteSiteAnalysersByChromosome.put(chromosome, new ImmutableIntervalTree<>(entries));
        }
    }

    public void processRead(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag())
            return;

        String chromosome = read.getReferenceName();
        int readAlignmentStart = read.getAlignmentStart();
        int readAlignmentEnd = read.getAlignmentEnd();

        ImmutableIntervalTree<MicrosatelliteSiteData> chromosomeAnalysers = mMicrosatelliteSiteAnalysersByChromosome.get(chromosome);
        if(chromosomeAnalysers == null || chromosomeAnalysers.isEmpty())
        {
            return;
        }

        // get analysers which are fully contained within the alignment of the read
        chromosomeAnalysers
                .containedIntervals(readAlignmentStart, readAlignmentEnd)
                .stream()
                .map(Pair::getValue)
                .forEach(analyser -> analyser.addReadToStats(read));
    }

    public Collection<MicrosatelliteSiteData> getMicrosatelliteSiteAnalysers()
    {
        return mMicrosatelliteSiteAnalysers;
    }
}
