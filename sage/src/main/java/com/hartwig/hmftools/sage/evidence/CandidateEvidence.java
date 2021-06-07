package com.hartwig.hmftools.sage.evidence;

import java.io.File;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.RefContextConsumer;
import com.hartwig.hmftools.sage.context.RefContextFactory;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.coverage.GeneCoverage;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class CandidateEvidence
{

    private static final Logger LOGGER = LogManager.getLogger(CandidateEvidence.class);

    private final SageConfig config;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panel;
    private final ReferenceSequenceFile refGenome;
    private final SamSlicerFactory samSlicerFactory;
    private final Coverage coverage;

    public CandidateEvidence(@NotNull final SageConfig config, @NotNull final List<VariantHotspot> hotspots, final List<GenomeRegion> panel,
            @NotNull final SamSlicerFactory samSlicerFactory, @NotNull final ReferenceSequenceFile refGenome, final Coverage coverage)
    {
        this.config = config;
        this.panel = panel;
        this.samSlicerFactory = samSlicerFactory;
        this.hotspots = hotspots;
        this.refGenome = refGenome;
        this.coverage = coverage;
    }

    @NotNull
    public List<AltContext> get(@NotNull final String sample, @NotNull final String bamFile, @NotNull final RefSequence refSequence,
            @NotNull final GenomeRegion bounds)
    {
        LOGGER.debug("Variant candidates {} position {}:{}", sample, bounds.chromosome(), bounds.start());
        final List<GeneCoverage> geneCoverage = coverage.coverage(sample, bounds.chromosome());
        final RefContextFactory candidates = new RefContextFactory(config, sample, hotspots, panel);
        final RefContextConsumer refContextConsumer = new RefContextConsumer(config, bounds, refSequence, candidates);

        final Consumer<SAMRecord> consumer = record ->
        {
            refContextConsumer.accept(record);
            if(!geneCoverage.isEmpty())
            {
                final GenomeRegion alignment =
                        GenomeRegions.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
                geneCoverage.forEach(x -> x.accept(alignment));
            }
        };

        return get(bamFile, bounds, consumer, candidates);
    }

    @NotNull
    private List<AltContext> get(@NotNull final String bamFile, @NotNull final GenomeRegion bounds,
            @NotNull final Consumer<SAMRecord> recordConsumer, @NotNull final RefContextFactory candidates)
    {
        final List<AltContext> altContexts = Lists.newArrayList();

        final SamSlicer slicer = samSlicerFactory.create(bounds);
        try(final SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(config.validationStringency())
                .referenceSource(new ReferenceSource(refGenome))
                .open(new File(bamFile)))
        {

            // First parse
            slicer.slice(tumorReader, recordConsumer);

            // Add all valid alt contexts
            altContexts.addAll(candidates.altContexts());
        } catch(Exception e)
        {
            throw new CompletionException(e);
        }

        return altContexts;
    }
}
