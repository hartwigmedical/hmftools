package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.amber.AmberGender.determineGender;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.WINDOW_SIZE;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.AmberQCFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.purple.segment.ExpectedBAF;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;

import org.apache.commons.cli.ParseException;

public class AmberData
{
    private static int DEFAULT_READ_DEPTH = 100;

    public final Multimap<Chromosome,AmberBAF> ChromosomeBafs;
    public final Multimap<Chromosome,PCFPosition> TumorSegments;
    public final Map<String,List<BaseRegion>> PcfRegions;

    public final int AverageTumorDepth;
    public final double Contamination;
    public final Gender PatientGender;

    private static final double MAX_SOMATIC_TOTAL_READ_COUNT_PROPORTION = 1.2;
    private static final double MIN_SOMATIC_TOTAL_READ_COUNT_PROPORTION = 0.8;

    public AmberData(
            final String sampleId, final String amberDirectory, final boolean germlineOnlyMode, final RefGenomeVersion refGenVersion)
            throws ParseException, IOException
    {
        String qcFile = AmberQCFile.generateFilename(amberDirectory, sampleId);
        if(!new File(qcFile).exists())
        {
            throw new ParseException("Unable to open Amber QC file: " + qcFile);
        }

        PPL_LOGGER.info("reading Amber QC from {}", qcFile);
        Contamination = AmberQCFile.read(qcFile).contamination();

        String amberFilename = AmberBAFFile.generateAmberFilenameForReading(amberDirectory, sampleId);
        if(!new File(amberFilename).exists())
        {
            throw new ParseException("Unable to open Amber BAF file: " + amberFilename);
        }

        PPL_LOGGER.info("reading Amber BAFs from {}", amberFilename);
        ChromosomeBafs = AmberBAFFile.read(amberFilename, !germlineOnlyMode);

        if(!germlineOnlyMode)
        {
            String pcfFilename = PCFFile.generateBAFFilename(amberDirectory, sampleId);

            if(!new File(pcfFilename).exists())
            {
                throw new ParseException("Unable to open Amber PCF file: " + pcfFilename);
            }

            PPL_LOGGER.info("reading Amber PCFs from {}", pcfFilename);

            TumorSegments = PCFFile.readPositions(WINDOW_SIZE, PCFSource.TUMOR_BAF, pcfFilename);
            PcfRegions = PCFFile.loadChrBaseRegions(pcfFilename);
        }
        else
        {
            TumorSegments = ArrayListMultimap.create();
            PcfRegions = Collections.emptyMap();
        }

        AverageTumorDepth = (int) Math.round(ChromosomeBafs.values()
                .stream()
                .mapToInt(AmberBAF::tumorDepth)
                .filter(x -> x > 0)
                .average()
                .orElse(DEFAULT_READ_DEPTH));

        PatientGender = determineGender(refGenVersion, ChromosomeBafs);

        PPL_LOGGER.info("Amber average tumor depth({}) ambiguous BAF({})",
                AverageTumorDepth, String.format("%.3f", ExpectedBAF.expectedBAF(AverageTumorDepth)));
    }

    public int minSomaticTotalReadCount()
    {
        return (int) Math.floor(MIN_SOMATIC_TOTAL_READ_COUNT_PROPORTION * AverageTumorDepth);
    }
    public int maxSomaticTotalReadCount()
    {
        return (int) Math.ceil(MAX_SOMATIC_TOTAL_READ_COUNT_PROPORTION * AverageTumorDepth);
    }

    public void clearCache()
    {
        ChromosomeBafs.clear();
        TumorSegments.clear();
    }

    @VisibleForTesting
    public AmberData(int averageTumorDepth, final Gender patientGender)
    {
        ChromosomeBafs = ArrayListMultimap.create();
        TumorSegments = ArrayListMultimap.create();
        PcfRegions = Maps.newHashMap();

        AverageTumorDepth = averageTumorDepth;
        Contamination = 0;
        PatientGender = patientGender;
    }
}
