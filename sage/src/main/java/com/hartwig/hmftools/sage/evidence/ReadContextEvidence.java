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
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.select.SamRecordSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class ReadContextEvidence {

    private static final Logger LOGGER = LogManager.getLogger(ReadContextEvidence.class);

    private final int minQuality;
    private final int typicalReadLength;
    private final SageConfig sageConfig;
    private final SamSlicerFactory samSlicerFactory;
    private final ReferenceSequenceFile refGenome;
    private final ReadContextCounterFactory factory;

    public ReadContextEvidence(@NotNull final SageConfig config, @NotNull final SamSlicerFactory samSlicerFactory,
            @NotNull final ReferenceSequenceFile refGenome, final Map<String, QualityRecalibrationMap> qualityRecalibrationMap) {
        this.minQuality = config.minMapQuality();
        this.sageConfig = config;
        this.samSlicerFactory = samSlicerFactory;
        this.refGenome = refGenome;
        this.factory = new ReadContextCounterFactory(config, qualityRecalibrationMap);
        this.typicalReadLength = config.typicalReadLength();
    }

    @NotNull
    public List<ReadContextCounter> get(@NotNull final List<Candidate> candidates, @NotNull final String sample,
            @NotNull final String bam) {

        final List<ReadContextCounter> counters = factory.create(sample, candidates);
        if (candidates.isEmpty()) {
            return counters;
        }

        final Candidate firstCandidate = candidates.get(0);
        final Candidate lastCandidate = candidates.get(candidates.size() - 1);

        final GenomeRegion bounds = GenomeRegions.create(firstCandidate.chromosome(),
                Math.max(firstCandidate.position() - typicalReadLength, 1),
                lastCandidate.position() + typicalReadLength);
        final SamSlicer slicer = samSlicerFactory.create(bounds);

        final SamRecordSelector<ReadContextCounter> consumerSelector = new SamRecordSelector<>(counters);

        try (final SamReader tumorReader = SamReaderFactory.makeDefault()
                .referenceSource(new ReferenceSource(refGenome))
                .open(new File(bam))) {
            slicer.slice(tumorReader, samRecord -> {

                if (samRecord.getMappingQuality() >= minQuality) {
                    int numberOfEvents = numberOfEvents(samRecord);
                    consumerSelector.select(samRecord, x -> x.accept(samRecord, sageConfig, numberOfEvents));
                }

            });
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return counters;
    }

    public int numberOfEvents(@NotNull final SAMRecord record) {
        Object nm = record.getAttribute("NM");
        if (!(nm instanceof Integer)) {
            return 0;
        }

        int additionalIndels = 0;
        for (CigarElement cigarElement : record.getCigar()) {
            switch (cigarElement.getOperator()) {
                case D:
                case I:
                    additionalIndels += (cigarElement.getLength() - 1);
            }
        }

        return (Integer) nm - additionalIndels;
    }

}
