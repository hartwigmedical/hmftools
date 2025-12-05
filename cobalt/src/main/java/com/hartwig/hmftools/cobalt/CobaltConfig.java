package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_GC_RATIO_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_GC_RATIO_MIN;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_PCF_GAMMA;
import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.common.bam.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.ExecutorService;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.count.BamReadCounter;
import com.hartwig.hmftools.cobalt.diploid.DiploidRegionLoader;
import com.hartwig.hmftools.cobalt.diploid.DiploidStatus;
import com.hartwig.hmftools.cobalt.exclusions.ExcludedRegionsFile;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.cobalt.targeted.TargetRegions;
import com.hartwig.hmftools.cobalt.targeted.TargetedRegionsNormalisationFile;
import com.hartwig.hmftools.cobalt.targeted.WholeGenome;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.cobalt.CobaltGcMedianFile;
import com.hartwig.hmftools.common.cobalt.CobaltMedianRatioFile;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class CobaltConfig
{
    public static final String TUMOR_ONLY_DIPLOID_BED = "tumor_only_diploid_bed";
    private static final String MIN_MAPPING_QUALITY = "min_quality";
    public static final String PCF_GAMMA = "pcf_gamma";

    public static final String TARGET_REGION_NORM_FILE = "target_region_norm_file";
    private static final String INCLUDE_DUPLICATES = "include_duplicates";

    public static final String GC_RATIO_MIN = "gc_ratio_min";
    public static final String GC_RATIO_MAX = "gc_ratio_max";

    private static final String SKIP_PCF_CALC = "skip_pcf_calc";
    public static final String USE_NEW_SEGMENTER = "use_new_segmenter";

    public final String ReferenceId;
    public final String ReferenceBamPath;

    public final String TumorId;
    public final String TumorBamPath;

    public final String RefGenomePath;
    public final RefGenomeVersion RefGenVersion;
    public final String GcProfilePath;

    public final int Threads;
    public final String OutputDir;

    public final int MinMappingQuality;
    public final int PcfGamma;

    public final ValidationStringency BamStringency;
    public final boolean IncludeDuplicates;
    public final boolean SkipPcfCalc;
    public final boolean UseNewSegmenter;

    public final String TumorOnlyDiploidBed;
    public final String TargetRegionNormFile;
    public List<ChrBaseRegion> mExcludedRegions;

    // debug
    public final SpecificRegions SpecificChrRegions;

    public static final Logger CB_LOGGER = LogManager.getLogger(CobaltConfig.class);

    public CobaltConfig(final ConfigBuilder configBuilder)
    {
        TumorId = configBuilder.getValue(TUMOR);
        TumorBamPath = configBuilder.getValue(TUMOR_BAM);

        ReferenceId = configBuilder.getValue(REFERENCE);
        ReferenceBamPath = configBuilder.getValue(REFERENCE_BAM);

        GcProfilePath = configBuilder.getValue(GC_PROFILE);

        TumorOnlyDiploidBed = configBuilder.getValue(TUMOR_ONLY_DIPLOID_BED);
        TargetRegionNormFile = configBuilder.getValue(TARGET_REGION_NORM_FILE);
        RefGenomePath = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);

        // set global constants
        CobaltConstants.GC_RATIO_MIN = configBuilder.getDecimal(GC_RATIO_MIN);
        CobaltConstants.GC_RATIO_MAX = configBuilder.getDecimal(GC_RATIO_MAX);

        MinMappingQuality = configBuilder.getInteger(MIN_MAPPING_QUALITY);
        PcfGamma = configBuilder.getInteger(PCF_GAMMA);
        IncludeDuplicates = configBuilder.hasFlag(INCLUDE_DUPLICATES);

        BamStringency = BamUtils.validationStringency(configBuilder);
        OutputDir = parseOutputDir(configBuilder);
        Threads = parseThreads(configBuilder);

        SkipPcfCalc = configBuilder.hasFlag(SKIP_PCF_CALC);
        UseNewSegmenter = configBuilder.hasFlag(USE_NEW_SEGMENTER);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        loadExcludedRegions();
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, TUMOR_DESC);
        configBuilder.addPath(TUMOR_BAM, false, TUMOR_BAM_DESC);

        configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC);
        configBuilder.addPath(REFERENCE_BAM, false, REFERENCE_BAM_DESC);

        configBuilder.addPath(GC_PROFILE, true, GC_PROFILE_DESC);
        configBuilder.addPath(REF_GENOME, false, REF_GENOME_CFG_DESC + ", required when using CRAM files");
        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC, null);

        registerCommonConfig(configBuilder);

        configBuilder.addPath(TUMOR_ONLY_DIPLOID_BED, false, "Diploid regions for tumor-only mode");
        configBuilder.addPath(TARGET_REGION_NORM_FILE, false, "Targeted regions normalisation file");

        configBuilder.addInteger(MIN_MAPPING_QUALITY, "Min map quality", DEFAULT_MIN_MAPPING_QUALITY);
        configBuilder.addInteger(PCF_GAMMA, "Gamma value for copy number PCF", DEFAULT_PCF_GAMMA);
        configBuilder.addFlag(INCLUDE_DUPLICATES, "Include duplicate reads in depth counts");
        configBuilder.addFlag(SKIP_PCF_CALC, "Skip final PCF output");
        configBuilder.addFlag(USE_NEW_SEGMENTER, "Use native segmenter");

        addSpecificChromosomesRegionsConfig(configBuilder);

        addOutputDir(configBuilder);
        addThreadOptions(configBuilder);
        addValidationStringencyOption(configBuilder);
        addLoggingOptions(configBuilder);
    }

    public static void registerCommonConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimal(GC_RATIO_MIN, "Restrict GC ratios to above minimum", DEFAULT_GC_RATIO_MIN);
        configBuilder.addDecimal(GC_RATIO_MAX, "Restrict GC ratios to below maximum", DEFAULT_GC_RATIO_MAX);
    }

    public void validate() throws Exception
    {
        if(ReferenceId == null)
        {
            if(ReferenceBamPath != null)
            {
                throw new Exception(String.format("%s option not allowed in tumor only mode", REFERENCE_BAM));
            }
        }
        else if(TumorOnlyDiploidBed != null)
        {
            throw new Exception(String.format("%s option is only allowed in tumor only mode", TUMOR_ONLY_DIPLOID_BED));
        }

        if(TumorId == null)
        {
            if(TumorBamPath != null)
            {
                throw new Exception(String.format("%s option not allowed in germline only mode", TUMOR_BAM));
            }
        }

        if(!checkCreateOutputDir(OutputDir))
        {
            throw new Exception(String.format("failed to create output directory(%s)", OutputDir));
        }

        if(GcProfilePath.endsWith(".gz"))
        {
            throw new Exception(String.format("invalid GC-profile file(%s), must be uncompressed", GcProfilePath));
        }

        if(mExcludedRegions == null)
        {
            throw new Exception("Excluded regions could not be loaded");
        }
    }

    public SamReaderFactory readerFactory()
    {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(BamStringency);

        if(RefGenomePath != null && !RefGenomePath.isEmpty())
        {
            return readerFactory.referenceSource(new ReferenceSource(new File(RefGenomePath)));
        }

        return readerFactory;
    }

    public CobaltScope scope()
    {
        if(TargetRegionNormFile == null)
        {
            return new WholeGenome();
        }
        else
        {
            TargetedRegionsNormalisationFile enrichmentFile = new TargetedRegionsNormalisationFile(TargetRegionNormFile);
            final RefGenomeCoordinates refGenomeCoordinates = RefGenomeCoordinates.refGenomeCoordinates(RefGenVersion);
            TargetRegions.ChromosomeData chromosomeData = chromosome -> refGenomeCoordinates.lengths().get(chromosome);
            return new TargetRegions(enrichmentFile.load(), chromosomeData);
        }
    }

    public ListMultimap<Chromosome, GCProfile> gcProfileData()
    {
        try
        {
            return GCProfileFactory.loadGCContent(WINDOW_SIZE, GcProfilePath);
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to load GC profile data: {}", e.toString());
            return ArrayListMultimap.create();
        }
    }

    public ListMultimap<Chromosome, DiploidStatus> diploidRegions()
    {
        if(TumorOnlyDiploidBed == null)
        {
            return ArrayListMultimap.create();
        }
        try
        {
            return new DiploidRegionLoader(TumorOnlyDiploidBed).regions();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to load diploid regions bed file: {}", e.toString());
            return ArrayListMultimap.create();
        }
    }

    public BamReadCounter tumorBamReader(ExecutorService executorService) throws IOException
    {
        return bamReader(executorService, TumorBamPath);
    }

    public BamReadCounter referenceBamReader(ExecutorService executorService) throws IOException
    {
        return bamReader(executorService, ReferenceBamPath);
    }

    public String cobaltRatiosFileName()
    {
        return CobaltRatioFile.generateFilename(OutputDir, TumorId != null ? TumorId : ReferenceId);
    }

    public String tumorPcfFileName()
    {
        return PCFFile.generateRatioFilename(OutputDir, TumorId);
    }

    public String referencePcfFileName()
    {
        return PCFFile.generateRatioFilename(OutputDir, ReferenceId);
    }

    public String medianRatiosFileName()
    {
        return CobaltMedianRatioFile.generateFilename(OutputDir, ReferenceId);
    }

    public String tumorGcMedianFileName()
    {
        return CobaltGcMedianFile.generateFilename(OutputDir, TumorId);
    }

    public String referenceGcMedianFileName()
    {
        return CobaltGcMedianFile.generateFilename(OutputDir, ReferenceId);
    }

    public List<ChrBaseRegion> excludedRegions()
    {
        return mExcludedRegions;
    }

    public RefGenomeVersion genome()
    {
        return RefGenVersion;
    }

    private void loadExcludedRegions()
    {
        try
        {
            String version = RefGenVersion.is37() ? "37" : "38";
            String refFile = "/regions/excluded_regions_v" + version + ".tsv";
            final InputStream resourceStream = Objects.requireNonNull(CobaltConfig.class.getResourceAsStream(refFile));
            final BufferedReader reader = new BufferedReader(new InputStreamReader(resourceStream));
            mExcludedRegions = new ExcludedRegionsFile(reader).regions();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to load excluded regions: {}", e.toString());
        }
    }

    private BamReadCounter bamReader(ExecutorService executorService, String bamPath) throws IOException
    {
        if(bamPath != null)
        {
            return new BamReadCounter(WINDOW_SIZE, this, executorService, bamPath);
        }
        return null;
    }
}
