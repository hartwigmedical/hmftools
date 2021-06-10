package com.hartwig.hmftools.svtools.fusion_likelihood;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.linx.LinxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.svtools.fusion_likelihood.CohortExpFusions.BUCKET_MAX;
import static com.hartwig.hmftools.svtools.fusion_likelihood.CohortExpFusions.BUCKET_MIN;
import static com.hartwig.hmftools.svtools.fusion_likelihood.CohortExpFusions.GENE_PAIR_DELIM;
import static com.hartwig.hmftools.svtools.fusion_likelihood.CohortExpFusions.GENOME_BASE_COUNT;
import static com.hartwig.hmftools.svtools.fusion_likelihood.CohortExpFusions.LONG_DDI_BUCKET;
import static com.hartwig.hmftools.svtools.fusion_likelihood.CohortExpFusions.MIN_BUCKET_LENGTH;
import static com.hartwig.hmftools.svtools.fusion_likelihood.CohortExpFusions.MIN_FUSION_RATE;
import static com.hartwig.hmftools.svtools.fusion_likelihood.CohortExpFusions.SHORT_INV_BUCKET;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseRegion.hasAnyPhaseMatch;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseRegion.regionsPhaseMatched;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_0;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_1;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_2;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_5P_UTR;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_MAX;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.typeAsInt;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_LONG_SAME_ARM;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_MEDIUM_INV;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_REMOTE;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_SHORT_INV;
import static com.hartwig.hmftools.svtools.fusion_likelihood.LikelihoodCalc.calcGeneOverlapAreas;
import static com.hartwig.hmftools.svtools.fusion_likelihood.LikelihoodCalc.reportGeneOverlaps;
import static com.hartwig.hmftools.svtools.fusion_likelihood.RegionAllocator.DEFAULT_BUCKET_REGION_RATIO;
import static com.hartwig.hmftools.svtools.fusion_likelihood.RegionAllocator.MIN_BLOCK_SIZE;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class FusionLikelihood
{
    private EnsemblDataCache mGeneTransCache;

    private final CohortExpFusions mCohortCalculator;

    // bucket demarcations - the actual buckets will be a set of consecutive pairs - eg length 0 -> length 1 etc
    private final List<Integer> mProximateBucketLengths;

    private final List<String> mRestrictedChromosomes;
    private final List<String> mRestrictedGeneIds;

    private final List<String[]> mGeneFusionPairs;
    private String mOutputDir;

    private static final String DEL_DUP_BUCKET_LENGTHS = "fl_del_dup_bucket_lengths";
    private static final String SHORT_INV_BUCKET_LENGTH = "fl_inv_bucket_length";

    private static final String GENE_PAIR_FILE = "gene_pair_file";

    // for testing
    private static final String LIMITED_GENE_IDS = "limited_gene_ids";
    private static final String LIMITED_CHROMOSOMES = "limited_chromosomes";

    public static final Logger FLC_LOGGER = LogManager.getLogger(FusionLikelihood.class);

    public FusionLikelihood()
    {
        mCohortCalculator = new CohortExpFusions();

        mProximateBucketLengths = Lists.newArrayList();
        mRestrictedChromosomes = Lists.newArrayList();
        mRestrictedGeneIds = Lists.newArrayList();
        mGeneFusionPairs = Lists.newArrayList();
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(DEL_DUP_BUCKET_LENGTHS, true, "Semi-colon separated DEL bucket lengths");
        options.addOption(SHORT_INV_BUCKET_LENGTH, true, "INV bucket length");
        options.addOption(LIMITED_GENE_IDS, true, "List of geneIds to test with");
        options.addOption(LIMITED_CHROMOSOMES, true, "List of chromosomes to test with");
        options.addOption(GENE_PAIR_FILE, true, "List of gene-pairs to calculate likelihood for");
    }

    public void initialise(final CommandLine cmdLineArgs, final EnsemblDataCache geneTransCache)
    {
        mGeneTransCache = geneTransCache;

        if(cmdLineArgs.hasOption(DEL_DUP_BUCKET_LENGTHS))
        {
            setBucketLengths(cmdLineArgs.getOptionValue(DEL_DUP_BUCKET_LENGTHS), mProximateBucketLengths);
        }

        if(cmdLineArgs.hasOption(LIMITED_CHROMOSOMES))
        {
            mRestrictedChromosomes.addAll(Arrays.stream(cmdLineArgs.getOptionValue(LIMITED_CHROMOSOMES)
                    .split(";"))
                    .collect(Collectors.toList()));
        }

        mOutputDir = parseOutputDir(cmdLineArgs);

        mCohortCalculator.initialiseLengths(mProximateBucketLengths, mRestrictedChromosomes);

        if(cmdLineArgs.hasOption(LIMITED_GENE_IDS))
        {
            mRestrictedGeneIds.addAll(Arrays.stream(cmdLineArgs.getOptionValue(LIMITED_GENE_IDS)
                    .split(";"))
                    .collect(Collectors.toList()));
        }
        else if(cmdLineArgs.hasOption(GENE_PAIR_FILE))
        {
            final String genePairFile = cmdLineArgs.getOptionValue(GENE_PAIR_FILE);

            FLC_LOGGER.info("calculating fusion likelihood for gene-pairs in file: {}", genePairFile);

            loadGenePairs(genePairFile);

            for(final String[] genePair : mGeneFusionPairs)
            {
                final String geneIdUp = genePair[0];
                final String geneIdDown = genePair[1];

                if (!mRestrictedGeneIds.contains(geneIdUp))
                    mRestrictedGeneIds.add(geneIdUp);

                if (!mRestrictedGeneIds.contains(geneIdDown))
                    mRestrictedGeneIds.add(geneIdDown);
            }
        }

        boolean limitedLoading = !mRestrictedGeneIds.isEmpty();

        if(!mGeneTransCache.load(limitedLoading))
        {
            FLC_LOGGER.error("Ensembl data cache load failed, exiting");
            return;
        }

        mGeneTransCache.createGeneIdDataMap();

        if(limitedLoading)
        {
            mGeneTransCache.loadTranscriptData(mRestrictedGeneIds);
        }

        if(mRestrictedGeneIds.size() <= 2)
            setLogVerbose(true);
    }

    public void run()
    {
        if(!mGeneFusionPairs.isEmpty())
        {
            FLC_LOGGER.info("calculating fusion likelihood for {} gene-pairs", mGeneFusionPairs.size());

            calculateSpecificFusionLikelihood();
        }
        else
        {
            FLC_LOGGER.info("generating genome-wide fusion likelihood data");

            generateGlobalExpectedFusionCounts();
            // fusionLikelihood.generateGlobalStats(outputDir);
        }
    }

    @VisibleForTesting
    public void initialise(final EnsemblDataCache geneTransCache, final List<Integer> delDepLengths)
    {
        mGeneTransCache = geneTransCache;
        mProximateBucketLengths.addAll(delDepLengths);
    }

    public void setLogVerbose(boolean toggle)
    {
        mCohortCalculator.setLogVerbose(toggle);
    }

    private void setBucketLengths(final String lengthData, List<Integer> bucketLengths)
    {
        if(lengthData.contains(";"))
        {
            Arrays.stream(lengthData.split(";")).forEach(x -> bucketLengths.add(Integer.parseInt(x)));
        }
        else if(lengthData.contains("-exp-"))
        {
            String[] startEnds = lengthData.split("-exp-");
            int startLength = Integer.parseInt(startEnds[0]);
            int endLength = Integer.parseInt(startEnds[1]);

            // add a bucket from the min to the first specified length
            bucketLengths.add(MIN_BUCKET_LENGTH);

            int bucketLength = startLength;
            while(bucketLength <= endLength)
            {
                bucketLengths.add(bucketLength);
                bucketLength *= 2;
            }
        }
    }

    private void writeGeneLikelihoodData()
    {
        FLC_LOGGER.info("writing output files");

        writeGeneData();
        writeProximateFusionData();
    }

    private void writeGeneData()
    {
        FLC_LOGGER.info("writing gene fusion data");

        try
        {
            String outputFilename = mOutputDir + "GFL_GENE_DATA.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("GeneId,GeneName,Chromosome,Arm,GeneStart,GeneEnd,Strand,ProteinCoding");
            writer.write(",FivePrimeUTR,Phase0,Phase1,Phase2,NonCoding");
            writer.write(",FivePrimeUTRPG,Phase0PG,Phase1PG,Phase2PG");
            writer.write(",ShortInvRateUp,ShortInvRateDown,MedInvRateUp,MedInvRateDown");
            writer.write(",LongDDIRateUp,LongDDIRateDown,RemoteRateUp,RemoteRateDown");
            writer.newLine();

            // adjustment factors to convert overlap base count into rates
            double remoteFusionFactor = 1.0 / (GENOME_BASE_COUNT * GENOME_BASE_COUNT);
            double sameArmFusionFactor = 1.0 / mCohortCalculator.getArmLengthFactor();
            int maxBucketLength = mCohortCalculator.getMaxBucketLength();
            double mediumInvFusionFactor = 1.0 / ((maxBucketLength - SHORT_INV_BUCKET) * GENOME_BASE_COUNT);
            double shortInvFusionFactor = 1.0 / (SHORT_INV_BUCKET * GENOME_BASE_COUNT);

            for(Map.Entry<String, List<GeneRangeData>> entry : mCohortCalculator.getChrGeneRangeDataMap().entrySet())
            {
                for(final GeneRangeData geneData :entry.getValue())
                {
                    int[] phaseCounts = new int[PHASE_MAX];
                    int[] phaseCountsPreGene = new int[PHASE_MAX];

                    geneData.getPhaseRegions().stream().forEach(x -> x.populateLengthCounts(phaseCounts, false));
                    geneData.getPhaseRegions().stream().forEach(x -> x.populateLengthCounts(phaseCountsPreGene, true));

                    writer.write(String.format("%s,%s,%s,%s,%d,%d,%d,%s",
                            geneData.GeneData.GeneId, geneData.GeneData.GeneName, geneData.GeneData.Chromosome, geneData.Arm,
                            geneData.GeneData.GeneStart, geneData.GeneData.GeneEnd, geneData.GeneData.Strand, geneData.hasProteinCoding()));

                    writer.write(String.format(",%d,%d,%d,%d,%d",
                            phaseCounts[typeAsInt(PHASE_5P_UTR)], phaseCounts[typeAsInt(PHASE_0)], phaseCounts[typeAsInt(PHASE_1)],
                            phaseCounts[typeAsInt(PHASE_2)], phaseCounts[typeAsInt(PHASE_NON_CODING)]));

                    writer.write(String.format(",%d,%d,%d,%d",
                            phaseCountsPreGene[typeAsInt(PHASE_5P_UTR)], phaseCountsPreGene[typeAsInt(PHASE_0)],
                            phaseCountsPreGene[typeAsInt(PHASE_1)], phaseCountsPreGene[typeAsInt(PHASE_2)]));

                    writer.write(String.format(",%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                            geneData.getBaseOverlapCountUpstream(NON_PROX_TYPE_SHORT_INV) * shortInvFusionFactor,
                            geneData.getBaseOverlapCountDownstream(NON_PROX_TYPE_SHORT_INV) * shortInvFusionFactor,
                            geneData.getBaseOverlapCountUpstream(NON_PROX_TYPE_MEDIUM_INV) * mediumInvFusionFactor,
                            geneData.getBaseOverlapCountDownstream(NON_PROX_TYPE_MEDIUM_INV) * mediumInvFusionFactor,
                            geneData.getBaseOverlapCountUpstream(NON_PROX_TYPE_LONG_SAME_ARM) * sameArmFusionFactor,
                            geneData.getBaseOverlapCountDownstream(NON_PROX_TYPE_LONG_SAME_ARM) * sameArmFusionFactor,
                            geneData.getBaseOverlapCountUpstream(NON_PROX_TYPE_REMOTE) * remoteFusionFactor,
                            geneData.getBaseOverlapCountDownstream(NON_PROX_TYPE_REMOTE) * remoteFusionFactor));

                    writer.newLine();
                }
            }

            closeBufferedWriter(writer);
        }
        catch (final IOException e)
        {
            FLC_LOGGER.error("error writing gene range data: {}", e.toString());
        }
    }

    private void writeProximateFusionData()
    {
        FLC_LOGGER.info("total gene-pair candidate count: dels({}) dups({})",
                mCohortCalculator.getDelGenePairCounts().size(), mCohortCalculator.getDupGenePairCounts().size());

        try
        {
            String outputFilename = mOutputDir + "GFL_DEL_DUP_PROXIMATES.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("Type,LengthMin,LengthMax,GeneIdUp,GeneNameUp,GeneIdDown,GeneNameDown,Chromosome,Strand,ProximateRate");
            writer.newLine();

            for(int i = 0; i <= 1; ++i)
            {
                boolean isDel = (i == 0);
                Map<String, Map<Integer,Long>> genePairCounts = isDel ? mCohortCalculator.getDelGenePairCounts() : mCohortCalculator.getDupGenePairCounts();

                for (Map.Entry<String,Map<Integer,Long>> entry : genePairCounts.entrySet())
                {
                    final String[] genePair = entry.getKey().split(GENE_PAIR_DELIM);
                    final String geneIdLower = genePair[0];
                    final String geneIdUpper = genePair[1];

                    EnsemblGeneData geneUp = null;
                    EnsemblGeneData geneDown = null;

                    Map<Integer,Long> bucketLengthCounts = entry.getValue();

                    for (Map.Entry<Integer,Long> bEntry : bucketLengthCounts.entrySet())
                    {
                        long overlapCount = bEntry.getValue();

                        int bucketIndex = bEntry.getKey();

                        int[] bucketMinMax = mCohortCalculator.getBucketLengthMinMax(bucketIndex);
                        int bucketWidth = bucketMinMax[BUCKET_MAX] - bucketMinMax[BUCKET_MIN];

                        double fusionRate = overlapCount / (bucketWidth * GENOME_BASE_COUNT);

                        if(fusionRate < MIN_FUSION_RATE)
                            continue;

                        if(geneUp == null && geneDown == null)
                        {
                            EnsemblGeneData geneLower = mGeneTransCache.getGeneDataById(geneIdLower);
                            EnsemblGeneData geneUpper = mGeneTransCache.getGeneDataById(geneIdUpper);
                            boolean isForwardStrand = (geneLower.Strand == 1);

                            geneUp = (isDel == isForwardStrand) ? geneLower : geneUpper;
                            geneDown = (!isDel == isForwardStrand) ? geneLower : geneUpper;
                        }

                        writer.write(String.format("%s,%d,%d,%s,%s,%s,%s,%s,%d,%.9f",
                                isDel ? "DEL" : "DUP", bucketMinMax[BUCKET_MIN], bucketMinMax[BUCKET_MAX],
                                geneUp.GeneId, geneUp.GeneName, geneDown.GeneId, geneDown.GeneName,
                                geneDown.Chromosome, geneDown.Strand, fusionRate));

                        writer.newLine();
                    }
                }
            }

            closeBufferedWriter(writer);
        }
        catch (final IOException e)
        {
            FLC_LOGGER.error("error writing gene-pair fusion candidates: {}", e.toString());
        }
    }

    private void generateGlobalExpectedFusionCounts()
    {
        mCohortCalculator.generateExpectedFusions(mGeneTransCache, mRestrictedChromosomes, mRestrictedGeneIds);
        writeGeneLikelihoodData();
    }

    public void generateGlobalStats()
    {
        mCohortCalculator.generateExpectedFusions(mGeneTransCache, mRestrictedChromosomes, mRestrictedGeneIds);
        reportGeneOverlaps(mCohortCalculator.getChrGeneRangeDataMap());
        mCohortCalculator.logGlobalCounts();
    }

    private void loadGenePairs(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
            return;

        // expected format: GeneIdUp,GeneIdDown
        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            // final String header = fileContents.get(0);
            fileContents.remove(0);

            fileContents.stream().map(x -> x.split(",")).filter(x -> x.length == 2).forEach(x -> mGeneFusionPairs.add(x));

            FLC_LOGGER.info("loaded {} specific gene fusions", mGeneFusionPairs.size());
        }
        catch (IOException e)
        {
            FLC_LOGGER.error("Failed to read specific gene fusions CSV file({}): {}", filename, e.toString());
        }
    }

    private void calculateSpecificFusionLikelihood()
    {
        // attempt to load cached gene phase data from file to avoid re-creating it each time

        // generate phase data for all required genes, which will also care of any same-gene fusion likelihoods
        final Map<String, GeneRangeData> geneIdRangeDataMap = mCohortCalculator.generateGeneRangeData(mGeneTransCache, mRestrictedGeneIds);

        try
        {
            String outputFilename = mOutputDir + "GFL_FUSION_PAIR_LIKELIHOOD.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("GeneIdUp,GeneNameUp,GeneIdDown,GeneNameDown,Type,LengthMin,LengthMax,FusionRate,GenePairRate");
            writer.newLine();

            int proximateLimit = mCohortCalculator.getMaxBucketLength();

            int bucketLengths = mProximateBucketLengths.size() - 1;
            final RegionAllocator[] regionAllocators = new RegionAllocator[bucketLengths];

            for(int i = 0; i < bucketLengths; ++i)
            {
                int blockSize = mProximateBucketLengths.get(i) / DEFAULT_BUCKET_REGION_RATIO;
                blockSize = max(blockSize, MIN_BLOCK_SIZE);
                regionAllocators[i] = new RegionAllocator(blockSize);
            }

            double sameArmFusionFactor = 1.0 / mCohortCalculator.getArmLengthFactor();
            int maxBucketLength = mCohortCalculator.getMaxBucketLength();
            double mediumInvFusionFactor = 1.0 / ((maxBucketLength - SHORT_INV_BUCKET) * GENOME_BASE_COUNT);
            double shortInvFusionFactor = 1.0 / (SHORT_INV_BUCKET * GENOME_BASE_COUNT);

            for(int i = 0; i < mGeneFusionPairs.size(); ++i)
            {
                final String[] genePair = mGeneFusionPairs.get(i);

                if(i > 0 && (i % 100) == 0)
                {
                    FLC_LOGGER.info("processed {} gene pairs", i);
                }

                final GeneRangeData geneUp = geneIdRangeDataMap.get(genePair[0]);
                final GeneRangeData geneDown = geneIdRangeDataMap.get(genePair[1]);

                geneUp.clearOverlapCounts();
                geneDown.clearOverlapCounts();

                boolean sameGene = geneUp.GeneData.GeneId.equals(geneDown.GeneData.GeneId);
                boolean sameStrand = geneUp.GeneData.Strand == geneDown.GeneData.Strand;
                boolean sameChromosome = geneUp.GeneData.Chromosome.equals(geneDown.GeneData.Chromosome);

                if(!sameGene)
                {
                    geneUp.setRestrictedStream(true);
                    geneDown.setRestrictedStream(false);
                }
                else
                {
                    geneUp.setRestrictedStream(null);
                    geneDown.setRestrictedStream(null);
                }

                final List<GeneRangeData> genePairList = Lists.newArrayList(geneUp, geneDown);

                if(sameChromosome)
                {
                    // work out whether these genes are proximate to each other or not
                    boolean proximateDelDup = false;

                    if(!sameGene && sameStrand)
                    {
                        int minDistance = min(abs(geneUp.GeneData.GeneStart - geneDown.GeneData.GeneStart),
                                abs(geneUp.GeneData.GeneEnd - geneDown.GeneData.GeneEnd));

                        minDistance = min(minDistance, abs(geneUp.GeneData.GeneStart - geneDown.GeneData.GeneEnd));
                        minDistance = min(minDistance, abs(geneUp.GeneData.GeneEnd - geneDown.GeneData.GeneStart));

                        if (minDistance <= proximateLimit)
                        {
                            proximateDelDup = true;
                        }
                    }

                    if(proximateDelDup || sameGene)
                    {
                        if (sameGene)
                        {
                            if(geneUp.hasProteinCoding())
                                mCohortCalculator.generateSameGeneCounts(geneUp);
                        }
                        else
                        {
                            mCohortCalculator.generateProximateCounts(genePairList, 1);
                            mCohortCalculator.generateProximateCounts(genePairList, -1);
                            mCohortCalculator.generateProximateCounts(genePairList, 0);
                        }

                        for(int j = 0; j <= 1; ++j)
                        {
                            boolean isDel = (j == 0);
                            Map<Integer,Long> proximateCounts = isDel ? geneUp.getDelFusionBaseCounts() : geneUp.getDupFusionBaseCounts();

                            for (Map.Entry<Integer,Long> bEntry : proximateCounts.entrySet())
                            {
                                long overlapArea = bEntry.getValue();

                                if(overlapArea == 0)
                                    continue;

                                int bucketIndex = bEntry.getKey();

                                int[] bucketMinMax = mCohortCalculator.getBucketLengthMinMax(bucketIndex);
                                int minBucketLen = bucketMinMax[BUCKET_MIN];
                                int maxBucketLen = bucketMinMax[BUCKET_MAX];
                                int bucketWidth = maxBucketLen - minBucketLen;

                                long geneOverlapArea = calcGeneOverlapAreas(
                                        geneUp, geneDown, isDel, minBucketLen, maxBucketLen, regionAllocators[bucketIndex]);

                                if(geneOverlapArea < overlapArea)
                                {
                                    // occurs rarely and not by more than 50%, not sure why
                                    double percent = overlapArea / (double)geneOverlapArea;
                                    FLC_LOGGER.debug("genes({} & {}) bucket({} -> {}) have inconsistent overlaps(fus={} gene={} perc={})",
                                            geneUp, geneDown, minBucketLen, maxBucketLen, overlapArea, geneOverlapArea, String.format("%.3f", percent));

                                    geneOverlapArea = overlapArea;
                                }

                                double genePairRate = geneOverlapArea / (bucketWidth * GENOME_BASE_COUNT);
                                double fusionRate = overlapArea / (bucketWidth * GENOME_BASE_COUNT);

                                writer.write(String.format("%s,%s,%s,%s,%s,%d,%d,%.6g,%.6g",
                                        geneUp.GeneData.GeneId, geneUp.GeneData.GeneName, geneDown.GeneData.GeneId, geneDown.GeneData.GeneName,
                                        isDel ? "DEL" : "DUP", bucketMinMax[BUCKET_MIN], bucketMinMax[BUCKET_MAX], fusionRate, genePairRate));

                                writer.newLine();
                            }
                        }
                    }
                    else
                    {
                        mCohortCalculator.generateNonProximateCounts(genePairList, 1);
                        mCohortCalculator.generateNonProximateCounts(genePairList, -1);
                        mCohortCalculator.generateNonProximateCounts(genePairList, 0);

                        if(geneUp.getBaseOverlapCountUpstream(NON_PROX_TYPE_SHORT_INV) > 0)
                        {

                            writer.write(String.format("%s,%s,%s,%s,%s,%d,%d,%.6g,0",
                                    geneUp.GeneData.GeneId, geneUp.GeneData.GeneName, geneDown.GeneData.GeneId, geneDown.GeneData.GeneName,
                                    "INV", MIN_BUCKET_LENGTH, SHORT_INV_BUCKET,
                                    geneUp.getBaseOverlapCountUpstream(NON_PROX_TYPE_SHORT_INV) * shortInvFusionFactor));

                            writer.newLine();
                        }

                        if(geneUp.getBaseOverlapCountUpstream(NON_PROX_TYPE_MEDIUM_INV) > 0)
                        {
                            writer.write(String.format("%s,%s,%s,%s,%s,%d,%d,%.6g,0",
                                    geneUp.GeneData.GeneId, geneUp.GeneData.GeneName, geneDown.GeneData.GeneId, geneDown.GeneData.GeneName,
                                    "INV", SHORT_INV_BUCKET, LONG_DDI_BUCKET,
                                    geneUp.getBaseOverlapCountUpstream(NON_PROX_TYPE_MEDIUM_INV) * mediumInvFusionFactor));

                            writer.newLine();
                        }

                        if(geneUp.getBaseOverlapCountUpstream(NON_PROX_TYPE_LONG_SAME_ARM) > 0)
                        {
                            double fusionRate = geneUp.getBaseOverlapCountUpstream(NON_PROX_TYPE_LONG_SAME_ARM) * sameArmFusionFactor;

                            writer.write(String.format("%s,%s,%s,%s,%s,%d,%d,%.6g,0",
                                    geneUp.GeneData.GeneId, geneUp.GeneData.GeneName, geneDown.GeneData.GeneId, geneDown.GeneData.GeneName,
                                    "LONG_DDI", LONG_DDI_BUCKET, LONG_DDI_BUCKET, fusionRate));

                            writer.newLine();
                        }
                    }
                }
                else
                {
                    // translocation fusion just require a check of overlapping bases regardless of any bucket length
                    int overlapArea = 0;

                    for (GenePhaseRegion regionUp : geneUp.getPhaseRegions())
                    {
                        for (GenePhaseRegion regionDown : geneDown.getPhaseRegions())
                        {
                            // test gene as an downstream vs all remote upstream phasing regions
                            if (hasAnyPhaseMatch(regionUp, regionDown, false) || regionsPhaseMatched(regionUp, regionDown))
                            {
                                overlapArea += regionUp.length() * regionDown.length();
                            }
                        }
                    }

                    if(overlapArea > 0)
                    {
                        double remoteFusionFactor = 1.0 / (GENOME_BASE_COUNT * GENOME_BASE_COUNT);

                        // GeneIdUp,GeneIdDown,Type,LengthMin,LengthMax,Likelihood
                        writer.write(String.format("%s,%s,%s,%s,%s,%d,%d,%.12f,0",
                                geneUp.GeneData.GeneId, geneUp.GeneData.GeneName, geneDown.GeneData.GeneId, geneDown.GeneData.GeneName,
                                "BND", 0, 0, overlapArea * remoteFusionFactor));

                        writer.newLine();
                    }
                }
            }

            closeBufferedWriter(writer);
        }
        catch (final IOException e)
        {
            FLC_LOGGER.error("error writing gene-pair fusion candidates: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addCmdLineArgs(options);
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log in verbose mode");
        options.addOption(ENSEMBL_DATA_DIR, true, "Ensembl gene transcript data cache directory");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if(cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        FusionLikelihood fusionLikelihood = new FusionLikelihood();

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(cmd.getOptionValue(ENSEMBL_DATA_DIR), RefGenomeVersion.V37);

        fusionLikelihood.initialise(cmd, ensemblDataCache);

        fusionLikelihood.run();

        FLC_LOGGER.info("gene-fusion likelihood calcs complete");
    }
}
