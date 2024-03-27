package com.hartwig.hmftools.common.basequal.jitter;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.collect.ImmutableIntervalTree;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMRecord;

public class SampleReadProcessor
{
    private final List<MicrosatelliteSiteAnalyser> mMicrosatelliteSiteAnalysers;
    private final Map<String, ImmutableIntervalTree<MicrosatelliteSiteAnalyser>> mMicrosatelliteSiteAnalysersByChromosome;

    public SampleReadProcessor(final Collection<RefGenomeMicrosatellite> refGenomeMicrosatellites)
    {
        mMicrosatelliteSiteAnalysers = refGenomeMicrosatellites.stream()
                .map(x -> new MicrosatelliteSiteAnalyser(x))
                .collect(Collectors.toList());

        mMicrosatelliteSiteAnalysersByChromosome = Maps.newHashMap();
        Multimap<String, MicrosatelliteSiteAnalyser> analysersByChromosome =
                Multimaps.index(mMicrosatelliteSiteAnalysers, analyser -> analyser.refGenomeMicrosatellite.chromosome());
        for(Map.Entry<String, Collection<MicrosatelliteSiteAnalyser>> chromosomeAnalysers : analysersByChromosome.asMap().entrySet())
        {
            String chromosome = chromosomeAnalysers.getKey();
            Collection<MicrosatelliteSiteAnalyser> analysers = chromosomeAnalysers.getValue();
            Collection<Pair<BaseRegion, MicrosatelliteSiteAnalyser>> entries = analysers
                    .stream()
                    .map(analyser -> Pair.of(analyser.refGenomeMicrosatellite.genomeRegion.baseRegion(), analyser))
                    .collect(Collectors.toList());

            mMicrosatelliteSiteAnalysersByChromosome.put(chromosome, new ImmutableIntervalTree<>(entries));
        }
    }

    public void processRead(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag())
        {
            return;
        }

        String chromosome = read.getReferenceName();
        int readAlignmentStart = read.getAlignmentStart();
        int readAlignmentEnd = read.getAlignmentEnd();

        ImmutableIntervalTree<MicrosatelliteSiteAnalyser> chromosomeAnalysers = mMicrosatelliteSiteAnalysersByChromosome.get(chromosome);
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

    public void processRead(final SAMRecord read, final ChrBaseRegion region)
    {
        if(read.getReadUnmappedFlag())
        {
            return;
        }

        if(!region.containsPosition(read.getAlignmentStart()))
        {
            return;
        }

        processRead(read);
    }

    public Collection<MicrosatelliteSiteAnalyser> getMicrosatelliteSiteAnalysers()
    {
        return mMicrosatelliteSiteAnalysers;
    }
}
