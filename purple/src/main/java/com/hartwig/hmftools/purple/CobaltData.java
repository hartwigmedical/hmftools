package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
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
    public final Map<Chromosome,List<PCFPosition>> TumorSegments;
    public final Map<Chromosome,List<PCFPosition>> ReferenceSegments;

    public final Gender gender()
    {
        return CobaltChromosomes.gender();
    }

    public CobaltData(
            final String referenceId, final String tumorId, final String cobaltDirectory,
            final Gender amberGender, final boolean tumorOnlyMode, final boolean germlineOnlyMode)
            throws ParseException, IOException
    {
        TumorSegments = Maps.newHashMap();
        ReferenceSegments = Maps.newHashMap();

        if(tumorId != null)
        {
            String pcfFileName = PCFFile.generateRatioFilename(cobaltDirectory, tumorId);
            PPL_LOGGER.info("reading Cobalt segments from {}", pcfFileName);
            TumorSegments.putAll(PCFFile.readPositions(PurpleConstants.WINDOW_SIZE, PCFSource.TUMOR_RATIO, pcfFileName));
        }

        if(referenceId != null)
        {
            String pcfFileName = PCFFile.generateRatioFilename(cobaltDirectory, referenceId);
            PPL_LOGGER.info("reading Cobalt segments from {}", pcfFileName);
            ReferenceSegments.putAll(PCFFile.readPositions(PurpleConstants.WINDOW_SIZE, PCFSource.REFERENCE_RATIO, pcfFileName));
        }

        String cobaltFilename = CobaltRatioFile.generateFilenameForReading(cobaltDirectory, tumorId);
        if(!new File(cobaltFilename).exists())
        {
            throw new ParseException("unable to open Cobalt ratio file: " + cobaltFilename);
        }

        PPL_LOGGER.info("reading Cobalt ratios from {}", cobaltFilename);
        Ratios = CobaltRatioFile.readWithGender(cobaltFilename, tumorOnlyMode ? amberGender : null, !germlineOnlyMode);

        List<MedianRatio> medianRatios = MedianRatioFactory.create(Ratios);
        CobaltChromosomes = new CobaltChromosomes(medianRatios, !tumorOnlyMode);
    }

    @VisibleForTesting
    public CobaltData(final CobaltChromosomes cobaltChromosomes)
    {
        CobaltChromosomes = cobaltChromosomes;
        Ratios = Maps.newHashMap();
        TumorSegments = Maps.newHashMap();
        ReferenceSegments = Maps.newHashMap();
    }
}
