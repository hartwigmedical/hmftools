package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.addGcProfilePath;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.hla.HlaCommon.hlaChromosome;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED_DESC;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.SampleDataFiles.GERMLINE_VARIANTS;
import static com.hartwig.hmftools.purple.germline.GermlineDeletionFrequency.COHORT_DEL_FREQ_FILE;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;
import com.hartwig.hmftools.purple.germline.GermlineDeletionFrequency;
import com.hartwig.hmftools.purple.region.ObservedRegionFactory;

import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ReferenceData
{
    public final RefGenomeVersion RefGenVersion;
    public final IndexedFastaSequenceFile RefGenome;

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

    private static final String TARGET_REGIONS_RATIOS = "target_regions_ratios";
    private static final String TARGET_REGION_MSI_INDELS = "target_regions_msi_indels";

    public ReferenceData(final ConfigBuilder configBuilder, final PurpleConfig config)
    {
        mIsValid = true;

        if(!configBuilder.hasValue(REF_GENOME) && !config.DriversOnly)
        {
            mIsValid = false;
            PPL_LOGGER.error(REF_GENOME + " is a mandatory argument");
        }

        final String refGenomePath = configBuilder.getValue(REF_GENOME);
        GcProfileFilename = configBuilder.getValue(GC_PROFILE);

        IndexedFastaSequenceFile refGenome = null;

        if(!config.DriversOnly)
        {
            try
            {
                refGenome = new IndexedFastaSequenceFile(new File(refGenomePath));
            }
            catch(Exception e)
            {
                mIsValid = false;
                PPL_LOGGER.error("failed to load ref genome: {}", e.toString());
            }
        }

        RefGenome = refGenome;

        RefGenVersion = RefGenomeVersion.from(configBuilder);
        PPL_LOGGER.info("using ref genome: {}", RefGenVersion);

        ChromosomeLengths = Maps.newHashMap();
        Centromeres = Maps.newHashMap();
        setChromosomeCoords();

        ObservedRegionFactory.setSpecificRegions(RefGenVersion, ChromosomeLengths, Centromeres);

        String somaticHotspotVcf = configBuilder.getValue(SOMATIC_HOTSPOT);
        String germlineHotspotVcf = configBuilder.getValue(GERMLINE_HOTSPOT);

        if(config.RunDrivers || config.DriversOnly)
        {
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
                    PPL_LOGGER.error("unable to open " + SOMATIC_HOTSPOT + " file " + somaticHotspotVcf);
                }
            }

            final List<DriverGene> driverGenes = Lists.newArrayList();

            try
            {
                driverGenes.addAll(DriverGenePanelConfig.driverGenes(configBuilder));
            }
            catch(IOException e)
            {
                mIsValid = false;
                PPL_LOGGER.error("unable to load driver genes: {}", e.toString());
            }

            DriverGenes = DriverGenePanelFactory.create(driverGenes);

            if(configBuilder.hasValue(GERMLINE_VARIANTS))
            {
                if(germlineHotspotVcf.isEmpty())
                {
                    mIsValid = false;
                    PPL_LOGGER.error(GERMLINE_HOTSPOT + " is a mandatory argument when running drivers");
                }

                if(!new File(germlineHotspotVcf).exists())
                {
                    mIsValid = false;
                    PPL_LOGGER.error("unable to open " + GERMLINE_HOTSPOT + " file " + germlineHotspotVcf);
                }
            }
        }
        else
        {
            DriverGenes = DriverGenePanelFactory.empty();
        }

        OtherReportableTranscripts = Maps.newHashMap();
        GeneTransCache = new EnsemblDataCache(configBuilder);
        loadGeneTransCache();

        if(mIsValid && config.tumorOnlyMode())
            HlaCommon.populateGeneData(GeneTransCache.getChrGeneDataMap().get(hlaChromosome(RefGenVersion)));

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

        CohortGermlineDeletions = new GermlineDeletionFrequency(configBuilder.getValue(COHORT_DEL_FREQ_FILE));

        TargetRegions = new TargetRegionsData(
                configBuilder.getValue(TARGET_REGIONS_BED),
                configBuilder.getValue(TARGET_REGIONS_RATIOS),
                configBuilder.getValue(TARGET_REGION_MSI_INDELS));
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

            PPL_LOGGER.debug("loaded {} alternative transcripts from {} genes",
                    additionalTransNames.size(), OtherReportableTranscripts.size());

            GeneTransCache.setRequiredData(true, false, false, true);
            mIsValid &= GeneTransCache.load(true);
            mIsValid &= GeneTransCache.loadTranscriptData(Collections.emptyList(), additionalTransNames);
        }
        else
        {
            GeneTransCache.setRequiredData(true, false, false, true);
            mIsValid &= GeneTransCache.load(false);
        }
    }

    public boolean isValid() { return mIsValid && TargetRegions.isValid(); }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        addRefGenomeConfig(configBuilder, false);

        configBuilder.addConfigItem(SOMATIC_HOTSPOT, false, "Path to somatic hotspot VCF", "");
        configBuilder.addConfigItem(GERMLINE_HOTSPOT, false, "Path to germline hotspot VCF", "");
        addGcProfilePath(configBuilder, false);
        configBuilder.addPath(COHORT_DEL_FREQ_FILE, false, "Path to cohort germline deletions frequency file");
        configBuilder.addPath(TARGET_REGIONS_BED, false, TARGET_REGIONS_BED_DESC);
        configBuilder.addPath(TARGET_REGIONS_RATIOS, false, "Path to target regions ratios file");
        configBuilder.addPath(TARGET_REGION_MSI_INDELS, false, "Path to target regions MSI INDELs file");
        EnsemblDataCache.addEnsemblDir(configBuilder, true);
        DriverGenePanelConfig.addGenePanelOption(configBuilder, false);
    }

    private void setChromosomeCoords()
    {
        RefGenomeCoordinates coordinates = RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = RefGenVersion.versionedChromosome(chromosome.toString());
            ChromosomeLengths.put(chromosome, GenomePositions.create(chrStr, coordinates.Lengths.get(chromosome)));
            Centromeres.put(chromosome, GenomePositions.create(chrStr, coordinates.Centromeres.get(chromosome)));
        }
    }

    @VisibleForTesting
    public ReferenceData(final PurpleConfig config)
    {
        mIsValid = true;
        GcProfileFilename = null;
        RefGenome = null;
        RefGenVersion = V38;
        ChromosomeLengths = Maps.newHashMap();
        Centromeres = Maps.newHashMap();
        setChromosomeCoords();
        DriverGenes = DriverGenePanelFactory.empty();

        OtherReportableTranscripts = Maps.newHashMap();
        GeneTransCache = new EnsemblDataCache("", RefGenVersion);
        SomaticHotspots = ArrayListMultimap.create();
        GermlineHotspots = ArrayListMultimap.create();
        CohortGermlineDeletions = new GermlineDeletionFrequency(null);
        TargetRegions = new TargetRegionsData(null, null, null);
    }
}
