package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.ratio.DiploidRatioSupplier.calcDiploidRatioResults;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.diploid.DiploidRatioBuilder;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.cobalt.MedianRatioFile;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCountFile;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.bed.BEDFeature;

public class RatioSupplier
{
    private final String mTumorId;
    private final String mReferenceId;
    private final String mOutputDir;

    class GermlineRatioSupplier
    {
        // processing states
        private final GCRatioSupplier gcRatioSupplier;
        private final ListMultimap<Chromosome, ReadRatio> gcDiploidRatios;

        ListMultimap<Chromosome, ReadRatio> getGcRatios() { return gcRatioSupplier.gcRatios(); }

        GermlineRatioSupplier(final Multimap<Chromosome, GCProfile> gcProfiles,
                final Multimap<Chromosome, CobaltCount> readCounts) throws IOException
        {
            CB_LOGGER.info("Applying ratio gc normalization");
            gcRatioSupplier = new GCRatioSupplier(gcProfiles, readCounts, CobaltCount::referenceReadCount);

            assert(mReferenceId != null);

            final List<MedianRatio> medianRatios = MedianRatioFactory.createFromReadRatio(gcRatioSupplier.gcRatios());

            CB_LOGGER.info("Persisting {} gc read count and reference ratio medians to {}", mReferenceId, mOutputDir);
            final String referenceGCMedianFilename = GCMedianReadCountFile.generateFilename(mOutputDir, mReferenceId);
            final String ratioMedianFilename = MedianRatioFile.generateFilename(mOutputDir, mReferenceId);

            GCMedianReadCountFile.write(referenceGCMedianFilename, gcRatioSupplier.gcMedianReadCount());
            MedianRatioFile.write(ratioMedianFilename, medianRatios);

            final CobaltChromosomes chromosomes = new CobaltChromosomes(medianRatios);
            if (chromosomes.hasGermlineAberrations())
            {
                CB_LOGGER.info("Found evidence of germline chromosomal aberrations: " + chromosomes.germlineAberrations()
                        .stream()
                        .map(Enum::toString)
                        .collect(Collectors.joining(",")));
            }

            CB_LOGGER.info("Applying ratio diploid normalization");
            gcDiploidRatios = calcDiploidRatioResults(chromosomes, getGcRatios());
        }
    }

    class TumorRatioSupplier
    {
        private final GCRatioSupplier gcRatioSupplier;

        ListMultimap<Chromosome, ReadRatio> getGcRatios() { return gcRatioSupplier.gcRatios(); }

        TumorRatioSupplier(final Multimap<Chromosome, GCProfile> gcProfiles,
                final Multimap<Chromosome, CobaltCount> readCounts) throws IOException
        {
            CB_LOGGER.info("Applying ratio gc normalization");
            gcRatioSupplier = new GCRatioSupplier(gcProfiles, readCounts, CobaltCount::tumorReadCount);

            CB_LOGGER.info("Persisting {} gc read count to {}", mTumorId, mOutputDir);
            final String tumorGCMedianFilename = GCMedianReadCountFile.generateFilename(mOutputDir, mTumorId);
            GCMedianReadCountFile.write(tumorGCMedianFilename, gcRatioSupplier.gcMedianReadCount());
        }
    }

    public RatioSupplier(final String reference, final String tumor, final String outputDirectory)
    {
        mTumorId = tumor;
        mReferenceId = reference;
        mOutputDir = outputDirectory;
    }

    @NotNull
    public Multimap<Chromosome, CobaltRatio> tumorOnly(final List<BEDFeature> bedFile, final Multimap<Chromosome, GCProfile> gcProfiles,
            final Multimap<Chromosome, CobaltCount> readCounts) throws IOException
    {
        var tumorRatioSupplier = new TumorRatioSupplier(gcProfiles, readCounts);

        final ListMultimap<Chromosome, ReadRatio> diploidRatios = new DiploidRatioBuilder(bedFile).build();

        // merge this ratios together into one cobalt ratio
        Multimap<Chromosome, CobaltRatio> result = ArrayListMultimap.create();
        final GenomePositionSelector<ReadRatio> tumorGCRatioSelector = GenomePositionSelectorFactory.create(tumorRatioSupplier.getGcRatios());
        final GenomePositionSelector<ReadRatio> diploidGCRatioSelector = GenomePositionSelectorFactory.create(diploidRatios);

        for (Chromosome chromosome : readCounts.keySet())
        {
            for (CobaltCount count : readCounts.get(chromosome))
            {
                if (diploidGCRatioSelector.select(count).isPresent())
                {
                    final CobaltRatio ratio = ImmutableCobaltRatio.builder()
                            .from(count)
                            .tumorGCRatio(tumorGCRatioSelector.select(count).map(ReadRatio::ratio).orElse(-1D))
                            .referenceGCRatio(-1D)
                            .referenceGCDiploidRatio(-1D)
                            .build();

                    result.put(chromosome, ratio);
                }
            }
        }

        return result;
    }

    @NotNull
    public Multimap<Chromosome, CobaltRatio> germlineOnly(final Multimap<Chromosome, GCProfile> gcProfiles,
            final Multimap<Chromosome, CobaltCount> readCounts) throws IOException
    {
        var germlineRatioSupplier = new GermlineRatioSupplier(gcProfiles, readCounts);

        // merge this ratios together into one cobalt ratio
        Multimap<Chromosome, CobaltRatio> result = ArrayListMultimap.create();
        final GenomePositionSelector<ReadRatio> referenceGCRatioSelector = GenomePositionSelectorFactory.create(germlineRatioSupplier.getGcRatios());
        final GenomePositionSelector<ReadRatio> referenceGCDiploidRatioSelector =
                GenomePositionSelectorFactory.create(germlineRatioSupplier.gcDiploidRatios);

        for (Chromosome chromosome : readCounts.keySet())
        {
            for (CobaltCount count : readCounts.get(chromosome))
            {
                final CobaltRatio ratio = ImmutableCobaltRatio.builder()
                        .from(count)
                        .tumorGCRatio(-1D)
                        .referenceGCRatio(referenceGCRatioSelector.select(count).map(ReadRatio::ratio).orElse(-1D))
                        .referenceGCDiploidRatio(referenceGCDiploidRatioSelector.select(count).map(ReadRatio::ratio).orElse(-1D))
                        .build();

                result.put(chromosome, ratio);
            }
        }

        return result;
    }

    @NotNull
    public Multimap<Chromosome, CobaltRatio> tumorNormalPair(final Multimap<Chromosome, GCProfile> gcProfiles,
            final Multimap<Chromosome, CobaltCount> readCounts) throws IOException
    {
        var tumorRatioSupplier = new TumorRatioSupplier(gcProfiles, readCounts);
        var germlineRatioSupplier = new GermlineRatioSupplier(gcProfiles, readCounts);

        // merge this ratios together into one cobalt ratio
        Multimap<Chromosome, CobaltRatio> result = ArrayListMultimap.create();
        final GenomePositionSelector<ReadRatio> tumorGCRatioSelector = GenomePositionSelectorFactory.create(tumorRatioSupplier.getGcRatios());
        final GenomePositionSelector<ReadRatio> referenceGCRatioSelector = GenomePositionSelectorFactory.create(germlineRatioSupplier.getGcRatios());
        final GenomePositionSelector<ReadRatio> referenceGCDiploidRatioSelector =
                GenomePositionSelectorFactory.create(germlineRatioSupplier.gcDiploidRatios);

        for (Chromosome chromosome : readCounts.keySet())
        {
            for (CobaltCount count : readCounts.get(chromosome))
            {
                final CobaltRatio ratio = ImmutableCobaltRatio.builder()
                        .from(count)
                        .tumorGCRatio(tumorGCRatioSelector.select(count).map(ReadRatio::ratio).orElse(-1D))
                        .referenceGCRatio(referenceGCRatioSelector.select(count).map(ReadRatio::ratio).orElse(-1D))
                        .referenceGCDiploidRatio(referenceGCDiploidRatioSelector.select(count).map(ReadRatio::ratio).orElse(-1D))
                        .build();

                result.put(chromosome, ratio);
            }
        }

        return result;
    }
}