package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltUtils.toCommonChromosomeMap;
import static com.hartwig.hmftools.cobalt.ratio.DiploidRatioSupplier.calcDiploidRatioResults;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.count.ReadCount;
import com.hartwig.hmftools.cobalt.diploid.DiploidRatioLoader;
import com.hartwig.hmftools.cobalt.targeted.TargetRegionEnrichment;
import com.hartwig.hmftools.cobalt.targeted.TargetedRatioBuilder;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.cobalt.MedianRatioFile;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCountFile;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class RatioSupplier
{
    private final String mTumorId;
    private final String mReferenceId;
    private final String mOutputDir;

    private final Collection<Chromosome> mChromosomes;
    private final Multimap<Chromosome, GCProfile> mGcProfiles;
    @Nullable private final Multimap<Chromosome, ReadCount> mReferenceCounts;
    @Nullable private final Multimap<Chromosome, ReadCount> mTumorCounts;

    private TargetRegionEnrichment mTargetRegionEnrichment = null;

    static class SampleRatios
    {
        // processing states
        protected final GcNormalizedRatioBuilder gcNormalizedRatioBuilder;
        protected final RatioBuilder ratioBuilder;

        ArrayListMultimap<Chromosome, ReadRatio> getRatios() { return ratioBuilder.ratios(); }

        SampleRatios(
                final String sampleId,
                final Multimap<Chromosome, ReadCount> readCounts,
                final Multimap<Chromosome, GCProfile> gcProfiles,
                @Nullable TargetRegionEnrichment targetRegionEnrichment,
                final String outputDir) throws IOException
        {
            gcNormalizedRatioBuilder = new GcNormalizedRatioBuilder(gcProfiles, readCounts);

            if (targetRegionEnrichment != null)
            {
                ratioBuilder = new TargetedRatioBuilder(
                        targetRegionEnrichment.regions(), targetRegionEnrichment.regionEnrichment(), gcNormalizedRatioBuilder.ratios());
            }
            else
            {
                ratioBuilder = gcNormalizedRatioBuilder;
            }

            CB_LOGGER.info("Persisting {} gc read count to {}", sampleId, outputDir);
            final String tumorGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDir, sampleId);
            GCMedianReadCountFile.write(tumorGCMedianFilename, gcNormalizedRatioBuilder.gcMedianReadCount());
        }
    }

    static class GermlineRatios extends SampleRatios
    {
        // processing states
        private final ArrayListMultimap<Chromosome, ReadRatio> gcDiploidRatios;

        GermlineRatios(final String referenceId,
                final Multimap<Chromosome, ReadCount> readCounts,
                final Multimap<Chromosome, GCProfile> gcProfiles,
                @Nullable TargetRegionEnrichment targetRegionEnrichment,
                final Collection<Chromosome> chromosomes,
                final String outputDir) throws IOException
        {
            super(referenceId, readCounts, gcProfiles, targetRegionEnrichment, outputDir);

            // TODO: check this
            final List<MedianRatio> medianRatios = MedianRatioFactory.createFromReadRatio(toCommonChromosomeMap(getRatios()));

            CB_LOGGER.info("Persisting {} gc ratio medians to {}", referenceId, outputDir);
            final String ratioMedianFilename = MedianRatioFile.generateFilename(outputDir, referenceId);
            MedianRatioFile.write(ratioMedianFilename, medianRatios);

            CB_LOGGER.info("Applying ratio diploid normalization");
            gcDiploidRatios = calcDiploidRatioResults(chromosomes, getRatios(), medianRatios);
        }
    }

    public RatioSupplier(final String reference, final String tumor, final String outputDirectory,
            final Multimap<Chromosome, GCProfile> gcProfiles,
            final Collection<Chromosome> chromosomes,
            @Nullable final Multimap<Chromosome, ReadCount> referenceCounts,
            @Nullable final Multimap<Chromosome, ReadCount> tumorCounts)
    {
        mTumorId = tumor;
        mReferenceId = reference;
        mOutputDir = outputDirectory;
        mGcProfiles = gcProfiles;
        mChromosomes = chromosomes;
        mReferenceCounts = referenceCounts;
        mTumorCounts = tumorCounts;
    }
    
    public void setTargetRegionEnrichment(TargetRegionEnrichment targetRegionEnrichment)
    {
        mTargetRegionEnrichment = targetRegionEnrichment;
    }

    @NotNull
    public Multimap<Chromosome, CobaltRatio> tumorOnly(final String diploidBedFile) throws IOException
    {
        if (mTumorCounts == null)
        {
            CB_LOGGER.fatal("Tumor count should not be null");
            throw new RuntimeException("tumor count is null");
        }
        var tumorRatios = new SampleRatios(mTumorId, mTumorCounts, mGcProfiles, mTargetRegionEnrichment, mOutputDir);
        final ArrayListMultimap<Chromosome, ReadRatio> diploidRatios = new DiploidRatioLoader(mChromosomes, diploidBedFile).build();

        // merge this ratios together into one cobalt ratio
        return mergeRatios(
                ArrayListMultimap.create(), mTumorCounts,
                diploidRatios, tumorRatios.getRatios(), diploidRatios);
    }

    @NotNull
    public Multimap<Chromosome, CobaltRatio> germlineOnly() throws IOException
    {
        if (mReferenceCounts == null)
        {
            CB_LOGGER.fatal("Reference count should not be null");
            throw new RuntimeException("reference count is null");
        }
        var germlineRatios = new GermlineRatios(mReferenceId, mReferenceCounts, mGcProfiles, mTargetRegionEnrichment, mChromosomes, mOutputDir);
        return mergeRatios(
                mReferenceCounts, ArrayListMultimap.create(),
                germlineRatios.getRatios(), ArrayListMultimap.create(), germlineRatios.gcDiploidRatios);
    }

    @NotNull
    public Multimap<Chromosome, CobaltRatio> tumorNormalPair() throws IOException
    {
        if (mReferenceCounts == null)
        {
            CB_LOGGER.fatal("Reference count should not be null");
            throw new RuntimeException("reference count is null");
        }
        if (mTumorCounts == null)
        {
            CB_LOGGER.fatal("Tumor count should not be null");
            throw new RuntimeException("tumor count is null");
        }
        var tumorRatios = new SampleRatios(mTumorId, mTumorCounts, mGcProfiles, mTargetRegionEnrichment, mOutputDir);
        var germlineRatios = new GermlineRatios(mReferenceId, mReferenceCounts, mGcProfiles, mTargetRegionEnrichment, mChromosomes, mOutputDir);

        return mergeRatios(
                mReferenceCounts, mTumorCounts,
                germlineRatios.getRatios(), tumorRatios.getRatios(), germlineRatios.gcDiploidRatios);
    }

    // merge everything together
    @NotNull
    private static Multimap<Chromosome, CobaltRatio> mergeRatios(
            @NotNull final Multimap<Chromosome, ReadCount> referenceCounts,
            @NotNull final Multimap<Chromosome, ReadCount> tumorCounts,
            @NotNull final ArrayListMultimap<Chromosome, ReadRatio> referenceRatios,
            @NotNull final ArrayListMultimap<Chromosome, ReadRatio> tumorRatios,
            @NotNull final ArrayListMultimap<Chromosome, ReadRatio> referenceDiploidRatios)
    {
        final Multimap<Chromosome, CobaltRatio> result = ArrayListMultimap.create();

        // find all the chromosomes
        Set<Chromosome> chromosomes = Sets.newIdentityHashSet();

        chromosomes.addAll(referenceCounts.keySet());
        chromosomes.addAll(tumorCounts.keySet());
        chromosomes.addAll(referenceRatios.keySet());
        chromosomes.addAll(tumorRatios.keySet());
        chromosomes.addAll(referenceDiploidRatios.keySet());

        for(Chromosome chromosome : chromosomes)
        {
            // try to merge all 5 lists

            Collection<ReadCount> refCountList = referenceCounts.get(chromosome);
            Collection<ReadCount> tumorCountList = tumorCounts.get(chromosome);

            // filter out NaN ratios
            Collection<ReadRatio> referenceRatioList = referenceRatios.get(chromosome).stream()
                    .filter(readRatio -> !Double.isNaN(readRatio.ratio())).collect(
                    Collectors.toList());
            Collection<ReadRatio> tumorRatioList = tumorRatios.get(chromosome).stream()
                    .filter(readRatio -> !Double.isNaN(readRatio.ratio())).collect(
                            Collectors.toList());
            Collection<ReadRatio> diploidRatioList = referenceDiploidRatios.get(chromosome).stream()
                    .filter(readRatio -> !Double.isNaN(readRatio.ratio())).collect(
                            Collectors.toList());

            // get all positions and add to a map
            Map<Integer, ImmutableCobaltRatio.Builder> positionRatioBuilders = new HashMap<>();

            for (Collection<? extends GenomePosition> l : List.of(refCountList, tumorCountList, referenceRatioList, tumorRatioList, diploidRatioList))
            {
                for (GenomePosition genomePosition : l)
                {
                    // set all initial values to -1
                    positionRatioBuilders.computeIfAbsent(
                            genomePosition.position(),
                            k -> ImmutableCobaltRatio.builder().from(genomePosition)
                                        .referenceReadCount(-1)
                                        .tumorReadCount(-1)
                                        .referenceGCRatio(-1D)
                                        .tumorGCRatio(-1D)
                                        .referenceGCDiploidRatio(-1D));
                }
            }

            // populate the values
            refCountList.forEach(rc -> positionRatioBuilders.get(rc.position()).referenceReadCount(rc.readCount()));
            tumorCountList.forEach(rc -> positionRatioBuilders.get(rc.position()).tumorReadCount(rc.readCount()));
            referenceRatioList.forEach(readRatio -> positionRatioBuilders.get(readRatio.position()).referenceGCRatio(readRatio.ratio()));
            tumorRatioList.forEach(readRatio -> positionRatioBuilders.get(readRatio.position()).tumorGCRatio(readRatio.ratio()));
            diploidRatioList.forEach(readRatio -> positionRatioBuilders.get(readRatio.position()).referenceGCDiploidRatio(readRatio.ratio()));

            List<CobaltRatio> cobaltRatios = positionRatioBuilders.values().stream()
                    .map(ImmutableCobaltRatio.Builder::build)
                    .collect(Collectors.toList());
            result.putAll(chromosome, cobaltRatios);
        }
        return result;
    }
}