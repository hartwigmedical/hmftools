package com.hartwig.hmftools.linx.fusion_likelihood;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions.BUCKET_MAX;
import static com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions.BUCKET_MIN;
import static com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions.GENE_PAIR_DELIM;
import static com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions.GENOME_BASE_COUNT;
import static com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions.MIN_BUCKET_LENGTH;
import static com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions.MIN_FUSION_RATE;
import static com.hartwig.hmftools.linx.fusion_likelihood.CohortExpFusions.SHORT_INV_BUCKET;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_0;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_1;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_2;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_5P_UTR;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_MAX;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.typeAsInt;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_LONG_SAME_ARM;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_MEDIUM_INV;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_REMOTE;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_SHORT_INV;
import static com.hartwig.hmftools.linx.fusion_likelihood.LikelihoodCalc.calcNonProximateLikelihood;
import static com.hartwig.hmftools.linx.fusion_likelihood.LikelihoodCalc.reportGeneOverlaps;
import static com.hartwig.hmftools.linx.types.SvaConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.linx.types.SvaConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.linx.types.SvaConfig.LOG_DEBUG;
import static com.hartwig.hmftools.linx.types.SvaConfig.formOutputPath;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

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
    private SvGeneTranscriptCollection mGeneTransCache;

    private CohortExpFusions mCohortCalculator;

    // bucket demarcations - the actual buckets will be a set of consecutive pairs - eg length 0 -> length 1 etc
    private List<Long> mProximateBucketLengths;

    private List<String> mRestrictedChromosomes;
    private List<String> mRestrictedGeneIds;
    private boolean mLogVerbose;

    private static final String DEL_DUP_BUCKET_LENGTHS = "fl_del_dup_bucket_lengths";
    private static final String SHORT_INV_BUCKET_LENGTH = "fl_inv_bucket_length";
    private static final String LIMITED_GENE_IDS = "limited_gene_ids"; // for testing
    private static final String LIMITED_CHROMOSOMES = "limited_chromosomes"; // for testing
    private static final String USE_PHASE_CACHE = "use_phase_cache";

    private static final Logger LOGGER = LogManager.getLogger(FusionLikelihood.class);

    public FusionLikelihood()
    {
        mCohortCalculator = new CohortExpFusions();

        mProximateBucketLengths = Lists.newArrayList();
        mRestrictedChromosomes = Lists.newArrayList();
        mRestrictedGeneIds = Lists.newArrayList();

        mLogVerbose = false;
    }


    public static void addCmdLineArgs(Options options)
    {
        options.addOption(DEL_DUP_BUCKET_LENGTHS, true, "Semi-colon separated DEL bucket lengths");
        options.addOption(SHORT_INV_BUCKET_LENGTH, true, "INV bucket length");
        options.addOption(LIMITED_GENE_IDS, true, "List of geneIds to test with");
        options.addOption(LIMITED_CHROMOSOMES, true, "List of chromosomes to test with");
    }

    public void initialise(final CommandLine cmdLineArgs, final SvGeneTranscriptCollection geneTransCache)
    {
        mGeneTransCache = geneTransCache;

        if(cmdLineArgs.hasOption(DEL_DUP_BUCKET_LENGTHS))
        {
            setBucketLengths(cmdLineArgs.getOptionValue(DEL_DUP_BUCKET_LENGTHS), mProximateBucketLengths);
        }

        int shortInvBucket = SHORT_INV_BUCKET;

        if(cmdLineArgs.hasOption(SHORT_INV_BUCKET_LENGTH))
        {
            shortInvBucket = Integer.parseInt(cmdLineArgs.getOptionValue(SHORT_INV_BUCKET_LENGTH));
        }

        if(cmdLineArgs.hasOption(LIMITED_CHROMOSOMES))
        {
            mRestrictedChromosomes = Arrays.stream(cmdLineArgs.getOptionValue(LIMITED_CHROMOSOMES).split(";"))
                    .collect(Collectors.toList());
        }

        mCohortCalculator.initialise(mProximateBucketLengths, shortInvBucket);

    }

    public void setRestrictedGeneIds(final List<String> geneIds) { mRestrictedGeneIds.addAll(geneIds); }

    @VisibleForTesting
    public void initialise(final SvGeneTranscriptCollection geneTransCache, final List<Long> delDepLengths)
    {
        mGeneTransCache = geneTransCache;
        mProximateBucketLengths.addAll(delDepLengths);
    }

    public void setLogVerbose(boolean toggle)
    {
        mLogVerbose = toggle;
        mCohortCalculator.setLogVerbose(toggle);
    }

    private void setBucketLengths(final String lengthData, List<Long> bucketLengths)
    {
        if(lengthData.contains(";"))
        {
            Arrays.stream(lengthData.split(";")).forEach(x -> bucketLengths.add(Long.parseLong(x)));
        }
        else if(lengthData.contains("-exp-"))
        {
            String[] startEnds = lengthData.split("-exp-");
            long startLength = Long.parseLong(startEnds[0]);
            long endLength = Long.parseLong(startEnds[1]);

            // add a bucket from the min to the first specified length
            bucketLengths.add(MIN_BUCKET_LENGTH);

            long bucketLength = startLength;
            while(bucketLength <= endLength)
            {
                bucketLengths.add(bucketLength);
                bucketLength *= 2;
            }
        }
    }


    public void writeGeneLikelihoodData(final String outputDir)
    {
        LOGGER.info("writing output files");

        writeGeneData(outputDir);
        writeTopProximateFusionCandidates(outputDir);
    }

    private void writeGeneData(final String outputDir)
    {
        LOGGER.info("writing gene fusion data");

        try
        {
            String outputFilename = outputDir + "GFL_GENE_DATA.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("GeneId,GeneName,Chromosome,Arm,GeneStart,GeneEnd,Strand");
            writer.write(",FivePrimeUTR,Phase0,Phase1,Phase2,NonCoding");
            writer.write(",FivePrimeUTRPG,Phase0PG,Phase1PG,Phase2PG");
            writer.write(",ShortInvRateUp,ShortInvRateDown,MedInvRateUp,MedInvRateDown");
            writer.write(",LongDDIRateUp,LongDDIRateDown,RemoteRateUp,RemoteRateDown");
            writer.newLine();

            // adjustment factors to convert overlap base count into rates
            double remoteFusionFactor = 1.0 / (GENOME_BASE_COUNT * GENOME_BASE_COUNT);
            double sameArmFusionFactor = 1.0 / mCohortCalculator.getArmLengthFactor();
            long maxBucketLength = mCohortCalculator.getMaxBucketLength();
            double mediumInvFusionFactor = 1.0 / ((maxBucketLength - SHORT_INV_BUCKET) * GENOME_BASE_COUNT);
            double shortInvFusionFactor = 1.0 / (SHORT_INV_BUCKET * GENOME_BASE_COUNT);

            for(Map.Entry<String, List<GeneRangeData>> entry : mCohortCalculator.getChrGeneRangeDataMap().entrySet())
            {
                for(final GeneRangeData geneData :entry.getValue())
                {
                    long[] phaseCounts = new long[PHASE_MAX];
                    long[] phaseCountsPreGene = new long[PHASE_MAX];

                    geneData.getPhaseRegions().stream().forEach(x -> x.populateLengthCounts(phaseCounts, false));
                    geneData.getPhaseRegions().stream().forEach(x -> x.populateLengthCounts(phaseCountsPreGene, true));

                    writer.write(String.format("%s,%s,%s,%s,%d,%d,%d",
                            geneData.GeneData.GeneId, geneData.GeneData.GeneName, geneData.GeneData.Chromosome, geneData.Arm,
                            geneData.GeneData.GeneStart, geneData.GeneData.GeneEnd, geneData.GeneData.Strand));

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
            LOGGER.error("error writing gene range data: {}", e.toString());
        }
    }

    private void writeTopProximateFusionCandidates(final String outputDir)
    {
        LOGGER.info("total gene-pair candidate count: dels({}) dups({})",
                mCohortCalculator.getDelGenePairCounts().size(), mCohortCalculator.getDupGenePairCounts().size());

        try
        {
            String outputFilename = outputDir + "GFL_DEL_DUP_PROXIMATES.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("Type,LengthMin,LengthMax,GeneIdUp,GeneNameUp,GeneIdDown,GeneNameDown,Chromosome,Strand,ProximateRate");
            writer.newLine();

            for(int i = 0; i <= 1; ++i)
            {
                boolean isDel = (i == 0);
                Map<String, Map<Integer, Long>> genePairCounts = isDel ? mCohortCalculator.getDelGenePairCounts() : mCohortCalculator.getDupGenePairCounts();

                for (Map.Entry<String, Map<Integer, Long>> entry : genePairCounts.entrySet())
                {
                    final String genePair[] = entry.getKey().split(GENE_PAIR_DELIM);
                    final String geneIdLower = genePair[0];
                    final String geneIdUpper = genePair[1];

                    EnsemblGeneData geneUp = null;
                    EnsemblGeneData geneDown = null;

                    Map<Integer, Long> bucketLengthCounts = entry.getValue();

                    for (Map.Entry<Integer, Long> bEntry : bucketLengthCounts.entrySet())
                    {
                        long overlapCount = bEntry.getValue();

                        int bucketIndex = bEntry.getKey();

                        long[] bucketMinMax = mCohortCalculator.getBucketLengthMinMax(isDel, bucketIndex);
                        long bucketWidth = bucketMinMax[BUCKET_MAX] - bucketMinMax[BUCKET_MIN];

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
            LOGGER.error("error writing gene-pair fusion candidates: {}", e.toString());
        }
    }

    public void generateGlobalFusionCounts(final String outputDir)
    {
        mCohortCalculator.generateExpectedFusions(mGeneTransCache, mRestrictedChromosomes, mRestrictedGeneIds);
        // mCohortCalculator.generateGenePhasingCounts(mGeneTransCache, mRestrictedChromosomes, mRestrictedGeneIds);
        // reportGeneOverlaps(mCohortCalculator.getChrGeneRangeDataMap());
        //mCohortCalculator.generateProximateFusionCounts();
        //mCohortCalculator.generateNonProximateCounts();
        // mCohortCalculator.logGlobalCounts();
        writeGeneLikelihoodData(outputDir);
    }

    public void calculateSpecificFusionLikelihood(final String filename, final String outputDir)
    {
        // attempt to load cached gene phase data from file to avoid re-creating it each time
        if(!Files.exists(Paths.get(filename)))
            return;

        // expected format: GeneIdUp,GeneIdDown

        // List<GeneRangeData> geneDataList = Lists.newArrayList();

        List<String[]> geneFusionPairs = Lists.newArrayList();

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
                return;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if(items.length < 2)
                    continue;

                final String geneIdUp = items[0];
                final String geneIdDown = items[1];

                geneFusionPairs.add(items);

                if(!mRestrictedGeneIds.contains(geneIdUp))
                    mRestrictedGeneIds.add(geneIdUp);

                if(!mRestrictedGeneIds.contains(geneIdDown))
                    mRestrictedGeneIds.add(geneIdDown);

            }
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read specific gene fusions CSV file({}): {}", filename, e.toString());
            return;
        }

        LOGGER.info("calculating fusion-likelihood for {} gene-pairs", geneFusionPairs.size());

        // generate phase data for all required genes, which will also care of any same-gene fusion likelihoods
        mCohortCalculator.initialiseGeneIdRangeDataMap();
        // mCohortCalculator.generateGenePhasingCounts(mGeneTransCache, mRestrictedChromosomes, mRestrictedGeneIds);

        try
        {
            String outputFilename = outputDir + "GFL_FUSION_PAIR_LIKELIHOOD.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("GeneIdUp,GeneIdDown,Type,LengthMin,LengthMax,LIkelihood");
            writer.newLine();

            long proximateLimit = mCohortCalculator.getMaxBucketLength();

            for(final String[] genePair : geneFusionPairs)
            {
                final GeneRangeData geneUp = mCohortCalculator.findGeneRangeData(genePair[0]);
                final GeneRangeData geneDown = mCohortCalculator.findGeneRangeData(genePair[1]);

                // work out whether these genes are proximate to each other or not
                boolean areProximate = false;

                if(geneUp.GeneData.Chromosome.equals(geneDown.GeneData.Chromosome) && geneUp.GeneData.Strand == geneDown.GeneData.Strand)
                {
                    long minDistance = min(abs(geneUp.GeneData.GeneStart - geneDown.GeneData.GeneStart),
                            abs(geneUp.GeneData.GeneEnd - geneDown.GeneData.GeneEnd));

                    minDistance = min(minDistance, abs(geneUp.GeneData.GeneStart - geneDown.GeneData.GeneEnd));
                    minDistance = min(minDistance, abs(geneUp.GeneData.GeneEnd - geneDown.GeneData.GeneStart));

                    if(minDistance <= proximateLimit)
                    {
                        areProximate = true;
                    }
                }

                if(areProximate)
                {
                    /*
                    for (GenePhaseRegion region1 : geneUp.getCombinedPhaseRegions())
                    {
                        for (GenePhaseRegion region2 : geneDown.getCombinedPhaseRegions())
                        {
                            testProximatePhaseRegions(geneUp, geneDown, region1, region2);
                        }
                    }

                    writer.write(String.format("%s,%d,%d,%s,%s,%s,%s,%s,%d,%.9f",
                            isDel ? "DEL" : "DUP", bucketMinMax[BUCKET_MIN], bucketMinMax[BUCKET_MAX],
                            geneUp.GeneId, geneUp.GeneName, geneDown.GeneId, geneDown.GeneName,
                            geneDown.Chromosome, geneDown.Strand, fusionRate));

                    writer.newLine();
                    */

                }
                else
                {
                    /*
                    long overlapCount = calcNonProximateLikelihood(geneUp, geneDown);

                    writer.write(String.format("%s,%d,%d,%s,%s,%s,%s,%s,%d,%.9f",
                            isDel ? "DEL" : "DUP", bucketMinMax[BUCKET_MIN], bucketMinMax[BUCKET_MAX],
                            geneUp.GeneId, geneUp.GeneName, geneDown.GeneId, geneDown.GeneName,
                            geneDown.Chromosome, geneDown.Strand, fusionRate));

                    writer.newLine();
                    */

                }

            }

            closeBufferedWriter(writer);
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing gene-pair fusion candidates: {}", e.toString());
        }

    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addCmdLineArgs(options);
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log in verbose mode");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Ensembl gene transcript data cache directory");
        options.addOption(USE_PHASE_CACHE, false, "Create and use gene phase regions cache file");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if(cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String outputDir = formOutputPath(cmd.getOptionValue(DATA_OUTPUT_DIR));

        LOGGER.info("Generating gene likelihood data");

        FusionLikelihood fusionLikelihood = new FusionLikelihood();

        SvGeneTranscriptCollection ensemblDataCache = new SvGeneTranscriptCollection();
        ensemblDataCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));

        List<String> restrictedGeneIds = Lists.newArrayList();
        if(cmd.hasOption(LIMITED_GENE_IDS))
        {
            restrictedGeneIds = Arrays.stream(cmd.getOptionValue(LIMITED_GENE_IDS).split(";")).collect(Collectors.toList());
            fusionLikelihood.setRestrictedGeneIds(restrictedGeneIds);
        }

        boolean limitedLoading = !restrictedGeneIds.isEmpty();

        if(!ensemblDataCache.loadEnsemblData(limitedLoading))
        {
            LOGGER.error("Ensembl data cache load failed, exiting");
            return;
        }

        ensemblDataCache.createGeneIdDataMap();

        if(limitedLoading)
        {
            ensemblDataCache.loadEnsemblTranscriptData(restrictedGeneIds);
        }

        fusionLikelihood.initialise(cmd, ensemblDataCache);

        if(restrictedGeneIds.size() <= 2)
            fusionLikelihood.setLogVerbose(true);

        fusionLikelihood.generateGlobalFusionCounts(outputDir);

        LOGGER.info("Gene likelihood data generation complete");
    }
}
