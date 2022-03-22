package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.bed.NamedBedFile.readBedFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.SampleDataFiles.GERMLINE_VARIANTS;
import static com.hartwig.hmftools.purple.germline.GermlineDeletionFrequency.COHORT_DEL_FREQ_FILE;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLengthFactory;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.purple.germline.GermlineDeletionFrequency;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ReferenceData
{
    public final RefGenomeVersion RefGenVersion;
    public final IndexedFastaSequenceFile RefGenome;
    public final RefGenomeCoordinates RefGeCoordinates;

    public final Map<Chromosome,GenomePosition> ChromosomeLengths;

    public final Map<Chromosome,GenomePosition> Centromeres;

    public final EnsemblDataCache GeneTransCache;

    public final DriverGenePanel DriverGenes;
    public final Map<String,List<String>> OtherReportableTranscripts;
    public final GermlineDeletionFrequency CohortGermlineDeletions;

    public final ListMultimap<Chromosome,VariantHotspot> SomaticHotspots;
    public final ListMultimap<Chromosome,VariantHotspot> GermlineHotspots;

    public final String GcProfileFilename;
    public final TargetRegionsData TargetRegions;

    private boolean mIsValid;

    private static final String SOMATIC_HOTSPOT = "somatic_hotspots";
    private static final String GERMLINE_HOTSPOT = "germline_hotspots";
    private static final String GC_PROFILE = "gc_profile";
    private static final String TARGET_REGIONS_RATIOS = "target_regions_ratios";

    public ReferenceData(final CommandLine cmd, final PurpleConfig config)
    {
        mIsValid = true;

        if(!cmd.hasOption(REF_GENOME) && !config.DriversOnly)
        {
            mIsValid = false;
            PPL_LOGGER.error(REF_GENOME + " is a mandatory argument");
        }

        final String refGenomePath = cmd.getOptionValue(REF_GENOME);
        GcProfileFilename = cmd.getOptionValue(GC_PROFILE);

        // TO-DO is this really necessary
        final Map<Chromosome, GenomePosition> lengthPositions = Maps.newHashMap();

        IndexedFastaSequenceFile refGenome = null;

        if(!config.DriversOnly)
        {
            try
            {
                refGenome = new IndexedFastaSequenceFile(new File(refGenomePath));

                SAMSequenceDictionary sequenceDictionary = refGenome.getSequenceDictionary();
                if(sequenceDictionary == null)
                {
                    throw new ParseException("Supplied ref genome must have associated sequence dictionary");
                }

                lengthPositions.putAll(fromLengths(ChromosomeLengthFactory.create(refGenome.getSequenceDictionary())));
            }
            catch(Exception e)
            {
                mIsValid = false;
                PPL_LOGGER.error("failed to load ref genome: {}", e.toString());
            }
        }

        RefGenome = refGenome;

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;
        PPL_LOGGER.info("using ref genome: {}", RefGenVersion);

        RefGeCoordinates = RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        final Map<Chromosome, String> chromosomeNames =
                lengthPositions.entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey, x -> x.getValue().chromosome()));

        ChromosomeLengths = toPosition(RefGeCoordinates.Lengths, chromosomeNames);
        Centromeres = toPosition(RefGeCoordinates.Centromeres, chromosomeNames);

        String somaticHotspotVcf = cmd.getOptionValue(SOMATIC_HOTSPOT, Strings.EMPTY);
        String germlineHotspotVcf = cmd.getOptionValue(GERMLINE_HOTSPOT, Strings.EMPTY);

        if(config.RunDrivers || config.DriversOnly)
        {
            if(!DriverGenePanelConfig.isConfigured(cmd))
            {
                mIsValid = false;
                PPL_LOGGER.error(DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION + " is a mandatory argument when running drivers");
            }

            if(!config.germlineMode())
            {
                if(somaticHotspotVcf.isEmpty())
                {
                    mIsValid = false;
                    PPL_LOGGER.error(SOMATIC_HOTSPOT + " is a mandatory argument when running drivers");
                }

                if(!new File(somaticHotspotVcf).exists())
                {
                    mIsValid = false;
                    PPL_LOGGER.error("Unable to open " + SOMATIC_HOTSPOT + " file " + somaticHotspotVcf);
                }
            }

            final List<DriverGene> driverGenes = Lists.newArrayList();

            try
            {
                driverGenes.addAll(DriverGenePanelConfig.driverGenes(cmd));
            }
            catch(IOException e)
            {
                mIsValid = false;
                PPL_LOGGER.error("Unable to load driver genes: {}", e.toString());
            }

            DriverGenes = DriverGenePanelFactory.create(driverGenes);

            if(cmd.hasOption(GERMLINE_VARIANTS))
            {
                if(germlineHotspotVcf.isEmpty())
                {
                    mIsValid = false;
                    PPL_LOGGER.error(GERMLINE_HOTSPOT + " is a mandatory argument when running drivers");
                }

                if(!new File(germlineHotspotVcf).exists())
                {
                    mIsValid = false;
                    PPL_LOGGER.error("Unable to open " + GERMLINE_HOTSPOT + " file " + germlineHotspotVcf);
                }
            }
        }
        else
        {
            DriverGenes = DriverGenePanelFactory.empty();
        }

        OtherReportableTranscripts = Maps.newHashMap();
        GeneTransCache = new EnsemblDataCache(cmd.getOptionValue(ENSEMBL_DATA_DIR, ""), RefGenVersion);
        loadGeneTransCache();

        SomaticHotspots = ArrayListMultimap.create();
        GermlineHotspots = ArrayListMultimap.create();

        try
        {
            SomaticHotspots.putAll(somaticHotspotVcf.equals(Strings.EMPTY) ?
                    ArrayListMultimap.create() : VariantHotspotFile.readFromVCF(somaticHotspotVcf));

            GermlineHotspots.putAll(germlineHotspotVcf.equals(Strings.EMPTY) ?
                    ArrayListMultimap.create() : VariantHotspotFile.readFromVCF(germlineHotspotVcf));
        }
        catch (IOException e)
        {
            mIsValid = false;
            PPL_LOGGER.error("failed to load hotspots: {}", e.toString());
        }

        CohortGermlineDeletions = new GermlineDeletionFrequency(cmd.getOptionValue(COHORT_DEL_FREQ_FILE));

        TargetRegions = new TargetRegionsData(config.TargetRegionsBed, cmd.getOptionValue(TARGET_REGIONS_RATIOS));
    }

    private void loadGeneTransCache()
    {
        // load transcripts with any alts from the driver gene panel in mind
        for(DriverGene driverGene : DriverGenes.driverGenes())
        {
            if(!driverGene.additionalReportedTranscripts().isEmpty())
            {
                OtherReportableTranscripts.put(driverGene.gene(), driverGene.additionalReportedTranscripts());
            }
        }

        if(!OtherReportableTranscripts.isEmpty())
        {
            List<String> additionalTransNames = Lists.newArrayList();
            OtherReportableTranscripts.values().forEach(x -> additionalTransNames.addAll(x));

            PPL_LOGGER.info("loaded {} alternative transcripts from {} genes",
                    additionalTransNames.size(), OtherReportableTranscripts.size());

            GeneTransCache.setRequiredData(true, false, false, true);
            GeneTransCache.load(true);
            GeneTransCache.loadTranscriptData(Lists.newArrayList(), additionalTransNames);
        }
        else
        {
            GeneTransCache.setRequiredData(true, false, false, true);
            GeneTransCache.load(false);
        }
    }



    public boolean isValid() { return mIsValid; }

    public static void addOptions(final Options options)
    {
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);

        options.addOption(SOMATIC_HOTSPOT, true, "Path to somatic hotspot VCF");
        options.addOption(GERMLINE_HOTSPOT, true, "Path to germline hotspot VCF");
        options.addOption(GC_PROFILE, true, "Path to GC profile");
        options.addOption(COHORT_DEL_FREQ_FILE, true, "Path to cohort germline deletions frequency file");
        options.addOption(TARGET_REGIONS_RATIOS, true, "Path to target regions ratios file");
        EnsemblDataCache.addEnsemblDir(options);
        DriverGenePanelConfig.addGenePanelOption(false, options);
    }

    private static Map<Chromosome, GenomePosition> toPosition(final Map<Chromosome,Integer> longs, final Map<Chromosome, String> contigMap)
    {
        final Map<Chromosome, GenomePosition> result = Maps.newHashMap();

        for(Map.Entry<Chromosome, String> entry : contigMap.entrySet())
        {
            final Chromosome chromosome = entry.getKey();
            final String contig = entry.getValue();
            if(longs.containsKey(chromosome))
            {
                result.put(chromosome, GenomePositions.create(contig, longs.get(chromosome)));
            }
        }

        return result;
    }

    private static Map<Chromosome, GenomePosition> fromLengths(final Collection<ChromosomeLength> lengths)
    {
        return lengths.stream()
                .filter(x -> HumanChromosome.contains(x.chromosome()))
                .collect(Collectors.toMap(x -> HumanChromosome.fromString(x.chromosome()),
                        item -> GenomePositions.create(item.chromosome(), item.length())));
    }


}
