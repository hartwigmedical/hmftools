package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.ratio.DiploidRatioSupplier.calcDiploidRatioResults;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.diploid.DiploidRatioBuilder;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFactory;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.cobalt.MedianRatioFile;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCountFile;
import com.hartwig.hmftools.common.genome.gc.GCProfile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.bed.BEDFeature;

public class RatioSupplier
{
    private final String tumor;
    private final String reference;
    private final String outputDirectory;

    public RatioSupplier(final String reference, final String tumor, final String outputDirectory)
    {
        this.tumor = tumor;
        this.reference = reference;
        this.outputDirectory = outputDirectory;
    }

    @NotNull
    public Multimap<Chromosome, CobaltRatio> tumorOnly(final List<BEDFeature> bedFile,
            @NotNull final Multimap<Chromosome, GCProfile> gcProfiles, @NotNull final Multimap<Chromosome, CobaltCount> readCounts)
            throws IOException
    {
        CB_LOGGER.info("Applying ratio gc normalization");
        final GCRatioSupplier gcRatioSupplier = new GCRatioSupplier(gcProfiles, readCounts);
        final ListMultimap<Chromosome, ReadRatio> tumorGCRatio = gcRatioSupplier.tumorRatios();

        final ListMultimap<Chromosome, ReadRatio> diploidRatio = new DiploidRatioBuilder(bedFile).build();

        CB_LOGGER.info("Persisting gc read count to {}", outputDirectory);
        final String tumorGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDirectory, tumor);
        GCMedianReadCountFile.write(tumorGCMedianFilename, gcRatioSupplier.tumorGCMedianReadCount());

        return CobaltRatioFactory.merge(readCounts, diploidRatio, tumorGCRatio, diploidRatio);
    }

    @NotNull
    public Multimap<Chromosome, CobaltRatio> tumorNormalPair(@NotNull final Multimap<Chromosome, GCProfile> gcProfiles,
            @NotNull final Multimap<Chromosome, CobaltCount> readCounts) throws IOException
    {
        CB_LOGGER.info("Applying ratio gc normalization");
        final GCRatioSupplier gcRatioSupplier = new GCRatioSupplier(gcProfiles, readCounts);
        final ListMultimap<Chromosome, ReadRatio> tumorGCRatio = gcRatioSupplier.tumorRatios();
        final ListMultimap<Chromosome, ReadRatio> referenceGCRatio = gcRatioSupplier.referenceRatios();

        final List<MedianRatio> medianRatios = MedianRatioFactory.createFromReadRatio(referenceGCRatio);
        final CobaltChromosomes chromosomes = new CobaltChromosomes(medianRatios);
        if(chromosomes.hasGermlineAberrations())
        {
            CB_LOGGER.info("Found evidence of germline chromosomal aberrations: " + chromosomes.germlineAberrations()
                    .stream()
                    .map(Enum::toString)
                    .collect(Collectors.joining(",")));
        }

        CB_LOGGER.info("Applying ratio diploid normalization");
        final ListMultimap<Chromosome, ReadRatio> referenceGCDiploidRatio = calcDiploidRatioResults(chromosomes, referenceGCRatio);

        CB_LOGGER.info("Persisting gc read count and reference ratio medians to {}", outputDirectory);
        final String tumorGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDirectory, tumor);
        final String referenceGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDirectory, reference);
        final String ratioMedianFilename = MedianRatioFile.generateFilename(outputDirectory, reference);
        GCMedianReadCountFile.write(tumorGCMedianFilename, gcRatioSupplier.tumorGCMedianReadCount());
        GCMedianReadCountFile.write(referenceGCMedianFilename, gcRatioSupplier.referenceGCMedianReadCount());
        MedianRatioFile.write(ratioMedianFilename, medianRatios);

        return CobaltRatioFactory.merge(readCounts, referenceGCRatio, tumorGCRatio, referenceGCDiploidRatio);
    }
}

