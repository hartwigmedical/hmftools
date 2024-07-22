package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.WINDOW_SIZE;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Maps;
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

    public final Map<Chromosome,List<CobaltRatio>> Ratios;
    public final Multimap<Chromosome, PCFPosition> TumorSegments;
    public final Multimap<Chromosome, PCFPosition> ReferenceSegments;

    public final Gender gender() { return CobaltChromosomes.gender(); }

    public CobaltData(
            final String referenceId, final String tumorId, final String cobaltDirectory,
            final Gender amberGender, final boolean tumorOnlyMode, final boolean germlineOnlyMode)
            throws ParseException, IOException
    {
        if(tumorId != null)
        {
            final String tumorSegmentFile = PCFFile.generateRatioFilename(cobaltDirectory, tumorId);
            if(!new File(tumorSegmentFile).exists())
            {
                throw new ParseException("unable to open Cobalt tumor pcf file: " + tumorSegmentFile);
            }

            PPL_LOGGER.info("reading Cobalt tumor segments from {}", tumorSegmentFile);
            TumorSegments = PCFFile.readPositions(WINDOW_SIZE, PCFSource.TUMOR_RATIO, tumorSegmentFile);
        }
        else
        {
            TumorSegments = ArrayListMultimap.create();
        }

        final String cobaltFilename = CobaltRatioFile.generateFilenameForReading(cobaltDirectory, tumorId);
        if(!new File(cobaltFilename).exists())
        {
            throw new ParseException("unable to open Cobalt ratio file: " + cobaltFilename);
        }

        PPL_LOGGER.info("reading Cobalt ratios from {}", cobaltFilename);
        Ratios = CobaltRatioFile.readWithGender(cobaltFilename, tumorOnlyMode ? amberGender : null, !germlineOnlyMode);

        if(referenceId != null)
        {
            final String referenceSegmentFile = PCFFile.generateRatioFilename(cobaltDirectory, referenceId);
            if(!new File(referenceSegmentFile).exists())
            {
                throw new ParseException("unable to open Cobalt reference PCF file: " + referenceSegmentFile);
            }

            PPL_LOGGER.info("reading Cobalt reference segments from {}", referenceSegmentFile);
            ReferenceSegments = PCFFile.readPositions(WINDOW_SIZE, PCFSource.REFERENCE_RATIO, referenceSegmentFile);
        }
        else
        {
            ReferenceSegments = ArrayListMultimap.create();
        }

        final List<MedianRatio> medianRatios = MedianRatioFactory.create(Ratios);
        CobaltChromosomes = new CobaltChromosomes(medianRatios, !tumorOnlyMode);
    }

    public void clearCache()
    {
        Ratios.clear();
        TumorSegments.clear();
        ReferenceSegments.clear();
    }

    @VisibleForTesting
    public CobaltData(final CobaltChromosomes cobaltChromosomes)
    {
        CobaltChromosomes = cobaltChromosomes;
        Ratios = Maps.newHashMap();
        TumorSegments = ArrayListMultimap.create();
        ReferenceSegments = ArrayListMultimap.create();
    }
}
