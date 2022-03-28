package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.WINDOW_SIZE;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.purple.segment.ExpectedBAF;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;

import org.apache.commons.cli.ParseException;

public class AmberData
{
    private static int DEFAULT_READ_DEPTH = 100;

    public final Gender PatientGender;

    public final Multimap<Chromosome, AmberBAF> ChromosomeBafs;

    public final Multimap<Chromosome, PCFPosition> TumorSegments;

    public final int AverageTumorDepth;

    public double Contamination;

    private static final double MAX_SOMATIC_TOTAL_READ_COUNT_PROPORTION = 1.4;
    private static final double MIN_SOMATIC_TOTAL_READ_COUNT_PROPORTION = 0.6;

    public int minSomaticTotalReadCount()
    {
        return (int) Math.floor(MIN_SOMATIC_TOTAL_READ_COUNT_PROPORTION * AverageTumorDepth);
    }

    public int maxSomaticTotalReadCount()
    {
        return (int) Math.ceil(MAX_SOMATIC_TOTAL_READ_COUNT_PROPORTION * AverageTumorDepth);
    }

    public AmberData(final String sampleId, final String amberDirectory, final boolean germlineOnlyMode) throws ParseException, IOException
    {
        final String qcFile = AmberQCFile.generateFilename(amberDirectory, sampleId);
        if(!new File(qcFile).exists())
        {
            throw new ParseException("Unable to open Amber QC file: " + qcFile);
        }

        PPL_LOGGER.info("reading Amber QC from {}", qcFile);
        Contamination = AmberQCFile.read(qcFile).contamination();

        final String amberFilename = AmberBAFFile.generateAmberFilenameForReading(amberDirectory, sampleId);
        if(!new File(amberFilename).exists())
        {
            throw new ParseException("Unable to open Amber BAF file: " + amberFilename);
        }

        PPL_LOGGER.info("reading Amber BAFs from {}", amberFilename);
        ChromosomeBafs = AmberBAFFile.read(amberFilename, !germlineOnlyMode);

        if(!germlineOnlyMode)
        {
            final String pcfFilename = PCFFile.generateBAFFilename(amberDirectory, sampleId);
            if(!new File(pcfFilename).exists())
            {
                throw new ParseException("Unable to open Amber PCF file: " + pcfFilename);
            }

            PPL_LOGGER.info("reading Amber PCFs from {}", pcfFilename);

            TumorSegments = PCFFile.readPositions(WINDOW_SIZE, PCFSource.TUMOR_BAF, pcfFilename);
        }
        else
        {
            TumorSegments = ArrayListMultimap.create();
        }

        AverageTumorDepth = (int) Math.round(ChromosomeBafs.values()
                .stream()
                .mapToInt(AmberBAF::tumorDepth)
                .filter(x -> x > 0)
                .average()
                .orElse(DEFAULT_READ_DEPTH));

        PPL_LOGGER.info("average Amber tumor depth is {} reads implying an ambiguous BAF of {}",
                AverageTumorDepth, String.format("%.3f", ExpectedBAF.expectedBAF(AverageTumorDepth)));

        PatientGender = Gender.fromAmber(ChromosomeBafs);
    }
}
