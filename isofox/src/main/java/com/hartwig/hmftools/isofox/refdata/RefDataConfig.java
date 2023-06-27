package com.hartwig.hmftools.isofox.refdata;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ER_FRAGMENT_LENGTHS;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LONG_FRAGMENT_LIMIT;
import static com.hartwig.hmftools.isofox.IsofoxConfig.READ_LENGTH;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_EXPECTED_RATE_LENGTHS;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MAX_FRAGMENT_SIZE;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesCommon.FL_LENGTH;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.IsofoxConstants;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.expression.ExpectedRatesCommon;

public class RefDataConfig
{
    public final int ReadLength;
    public final int MaxFragmentLength;
    public final List<FragmentSize> FragmentSizeData;
    public final boolean GenerateExpectedCounts;
    public final boolean GenerateGcRatios;
    public final RefGenomeVersion RefGenVersion;

    public final List<String> RestrictedGeneIds; // specific set of genes to process
    public final List<String> EnrichedGeneIds;

    public final int Threads;
    public final String OutputDir;
    public final String OutputId;

    private final static String GEN_EXPECTED_COUNTS = "expected_counts";
    private final static String GEN_GC_RATIOS = "expected_gc_ratios";

    public RefDataConfig(final ConfigBuilder configBuilder)
    {
        GenerateExpectedCounts = configBuilder.hasFlag(GEN_EXPECTED_COUNTS);
        GenerateGcRatios = configBuilder.hasFlag(GEN_GC_RATIOS);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        Threads = parseThreads(configBuilder);
        ReadLength = configBuilder.getInteger(READ_LENGTH);
        MaxFragmentLength = configBuilder.getInteger(LONG_FRAGMENT_LIMIT);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        EnrichedGeneIds = Lists.newArrayList();
        IsofoxConstants.populateEnrichedGeneIds(EnrichedGeneIds, RefGenVersion);

        FragmentSizeData = Lists.newArrayList();

        String[] fragLengths = configBuilder.getValue(ER_FRAGMENT_LENGTHS).split(ITEM_DELIM);
        for(int i = 0; i < fragLengths.length; ++i)
        {
            String[] flItem = fragLengths[i].split("-");
            int fragLength = Integer.parseInt(flItem[FL_LENGTH]);
            int fragFrequency = max(Integer.parseInt(flItem[ExpectedRatesCommon.FL_FREQUENCY]), 1);
            FragmentSizeData.add(new FragmentSize(fragLength, fragFrequency));
        }

        RestrictedGeneIds = Lists.newArrayList();
        if(configBuilder.hasValue(GENE_ID_FILE))
        {
            final String inputFile = configBuilder.getValue(GENE_ID_FILE);
            RestrictedGeneIds.addAll(loadGeneIdsFile(inputFile));

            if(!RestrictedGeneIds.isEmpty())
            {
                ISF_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
            }
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(GEN_EXPECTED_COUNTS, "Generate expected transcript counts");
        configBuilder.addFlag(GEN_GC_RATIOS, "Generate expected GC ratios");

        configBuilder.addInteger(LONG_FRAGMENT_LIMIT, "Max RNA fragment size", DEFAULT_MAX_FRAGMENT_SIZE);
        configBuilder.addRequiredInteger(READ_LENGTH, "Sample sequencing read length");

        configBuilder.addPath(GENE_ID_FILE, false, "Optional CSV file of genes to analyse");

        configBuilder.addConfigItem(
                ER_FRAGMENT_LENGTHS, false,
                "Fragment sizes and weights for expected transcript calcs (format: length1-freq1;length3-freq2 eg 100-10;150-20) in integer terms",
                DEFAULT_EXPECTED_RATE_LENGTHS);

        addRefGenomeConfig(configBuilder, false);

        addEnsemblDir(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }

    public RefDataConfig(final int readLength)
    {
        GenerateExpectedCounts = true;
        GenerateGcRatios = true;
        OutputDir = null;
        OutputId = null;
        Threads = 0;
        ReadLength = readLength;
        MaxFragmentLength = DEFAULT_MAX_FRAGMENT_SIZE;

        RefGenVersion = V37;

        EnrichedGeneIds = Lists.newArrayList();
        FragmentSizeData = Lists.newArrayList();
        RestrictedGeneIds = Lists.newArrayList();
    }
}
