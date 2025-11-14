package com.hartwig.hmftools.purple.data;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.io.IOException;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;
import com.hartwig.hmftools.purple.PurpleConstants;

import org.apache.commons.cli.ParseException;

public class PcfData
{
    private final File Directory;

    public PcfData(final String directory)
    {
        Directory = new File(directory);
        Preconditions.checkArgument(Directory.exists() , "Directory does not exist: " + directory);
        Preconditions.checkArgument(Directory.isDirectory(), "Not a directory: " + directory);
    }

    public ListMultimap<Chromosome, PCFPosition> loadCobaltSegments(String sample, PCFSource source) throws IOException, ParseException
    {
        if (sample == null)
        {
            return ArrayListMultimap.create();
        }
        String pcfFileName = PCFFile.generateRatioFilename(Directory.getAbsolutePath(), sample);
        if(!new File(pcfFileName).exists())
        {
            throw new org.apache.commons.cli.ParseException("unable to open Cobalt pcf file: " + pcfFileName);
        }

        PPL_LOGGER.info("reading Cobalt segments from {}", pcfFileName);
        return PCFFile.readPositions(PurpleConstants.WINDOW_SIZE, source, pcfFileName);
    }
}
