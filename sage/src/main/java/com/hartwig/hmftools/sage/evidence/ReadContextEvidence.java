package com.hartwig.hmftools.sage.evidence;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletionException;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.read.ReadContextCounterFactory;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.samtools.NumberEvents;
import com.hartwig.hmftools.sage.select.SamRecordSelector;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class ReadContextEvidence
{
    private final int mTypicalReadLength;
    private final SageConfig mSageConfig;
    private final ReferenceSequenceFile mRefGenome;
    private final ReadContextCounterFactory mFactory;

    public ReadContextEvidence(@NotNull final SageConfig config, @NotNull final ReferenceSequenceFile refGenome,
            final Map<String, QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mSageConfig = config;
        mRefGenome = refGenome;
        mFactory = new ReadContextCounterFactory(config, qualityRecalibrationMap);
        mTypicalReadLength = config.typicalReadLength();
    }

    @NotNull
    public List<ReadContextCounter> get(@NotNull final List<Candidate> candidates, @NotNull final String sample,
            @NotNull final String bam)
    {
        final List<ReadContextCounter> counters = mFactory.create(sample, candidates);
        if(candidates.isEmpty())
        {
            return counters;
        }

        final Candidate firstCandidate = candidates.get(0);
        final Candidate lastCandidate = candidates.get(candidates.size() - 1);

        final GenomeRegion bounds = GenomeRegions.create(firstCandidate.chromosome(),
                Math.max(firstCandidate.position() - mTypicalReadLength, 1),
                lastCandidate.position() + mTypicalReadLength);
        final SamSlicer slicer = new SamSlicer(0, bounds);

        final SamRecordSelector<ReadContextCounter> consumerSelector = new SamRecordSelector<>(counters);

        final RefSequence refSequence = new RefSequence(bounds, mRefGenome);

        try(final SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(mSageConfig.Stringency)
                .referenceSource(new ReferenceSource(mRefGenome))
                .open(new File(bam)))
        {
            slicer.slice(tumorReader, samRecord ->
            {
                int numberOfEvents = NumberEvents.numberOfEvents(samRecord, refSequence);
                consumerSelector.select(samRecord, x -> x.accept(samRecord, mSageConfig, numberOfEvents));

            });
        } catch(IOException e)
        {
            throw new CompletionException(e);
        }

        return counters;
    }
}
