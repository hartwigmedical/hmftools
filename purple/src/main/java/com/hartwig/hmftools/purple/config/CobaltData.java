package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.WINDOW_SIZE;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;

import org.apache.commons.cli.ParseException;

public class CobaltData
{
    public final CobaltChromosomes CobaltChromosomes;

    public final ListMultimap<Chromosome, CobaltRatio> Ratios;

    public final Multimap<Chromosome, PCFPosition> TumorSegments;

    public final Multimap<Chromosome, PCFPosition> ReferenceSegments;

    public final Gender gender() { return CobaltChromosomes.gender(); }

    public CobaltData(
            final String referenceId, final String tumorSample, final String cobaltDirectory,
            final Gender amberGender, final boolean tumorOnlyMode)
            throws ParseException, IOException
    {
        final String cobaltFilename = CobaltRatioFile.generateFilenameForReading(cobaltDirectory, tumorSample);
        if(!new File(cobaltFilename).exists())
        {
            throw new ParseException("uable to open cobalt ratio file: " + cobaltFilename);
        }

        String cobaltRefId = referenceId != null ? referenceId : CobaltRatioFile.TUMOR_ONLY_REFERENCE_SAMPLE;

        final String referenceSegmentFile = PCFFile.generateRatioFilename(cobaltDirectory, cobaltRefId);
        if(!new File(referenceSegmentFile).exists())
        {
            throw new ParseException("unable to open cobalt reference pcf file: " + referenceSegmentFile);
        }

        final String tumorSegmentFile = PCFFile.generateRatioFilename(cobaltDirectory, tumorSample);
        if(!new File(tumorSegmentFile).exists())
        {
            throw new ParseException("unable to open cobalt tumor pcf file: " + tumorSegmentFile);
        }

        PPL_LOGGER.info("reading cobalt ratios from {}", cobaltFilename);
        Ratios = tumorOnlyMode
                ? CobaltRatioFile.readTumorOnly(cobaltFilename, amberGender)
                : CobaltRatioFile.read(cobaltFilename);

        PPL_LOGGER.info("reading cobalt reference segments from {}", referenceSegmentFile);
        ReferenceSegments = PCFFile.readPositions(WINDOW_SIZE, PCFSource.REFERENCE_RATIO, referenceSegmentFile);

        PPL_LOGGER.info("reading cobalt tumor segments from {}", tumorSegmentFile);
        TumorSegments = PCFFile.readPositions(WINDOW_SIZE, PCFSource.TUMOR_RATIO, tumorSegmentFile);

        final List<MedianRatio> medianRatios = MedianRatioFactory.create(Ratios);
        CobaltChromosomes = new CobaltChromosomes(medianRatios);
    }
}
