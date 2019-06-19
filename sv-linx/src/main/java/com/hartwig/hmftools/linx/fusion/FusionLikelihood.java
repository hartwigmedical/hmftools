package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_LENGTHS;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArmLength;
import static com.hartwig.hmftools.linx.fusion.GenePhaseRegion.hasAnyPhaseMatch;
import static com.hartwig.hmftools.linx.fusion.GenePhaseRegion.mapExonPhase;
import static com.hartwig.hmftools.linx.fusion.GenePhaseRegion.regionsPhaseMatched;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_0;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_1;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_2;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_5P_UTR;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.NON_PROX_TYPE_LONG_SAME_ARM;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.NON_PROX_TYPE_REMOTE;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.NON_PROX_TYPE_SHORT_INV;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.PGD_DELIMITER;
import static com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection.nextTranscriptExons;
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
import com.google.common.collect.Maps;
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

    // bucket demarcations - the actual buckets will be a set of consecutive pairs - eg length 0 -> length 1 etc
    private List<Long> mDelBucketLengths;
    private List<Long> mDupBucketLengths;

    private final Map<String, List<GeneRangeData>> mChrForwardGeneDataMap;
    private final Map<String, List<GeneRangeData>> mChrReverseGeneDataMap;
    private Map<String,String> mGeneRangeDataMap; // populated from file

    private Map<String,Map<Integer,Long>> mDelGenePairCounts; // pair of gene-pairs to their overlap counts keyed by bucket index
    private Map<String,Map<Integer,Long>> mDupGenePairCounts;

    // cache the total counts (expressed as a length) for each unique phase combination per chromosome
    private Map<String, List<GenePhaseRegion>> mChrPhaseRegionsPosStrand;
    private Map<String, List<GenePhaseRegion>> mChrPhaseRegionsNegStrand;

    private List<String> mRestrictedChromosomes;
    private List<String> mRestrictedGeneIds;
    private long mArmLengthFactor;
    private boolean mLogVerbose;

    private static final String DEL_BUCKET_LENGTHS = "fl_del_bucket_lengths";
    private static final String DUP_BUCKET_LENGTHS = "fl_dup_bucket_lengths";
    private static final String SHORT_INV_BUCKET_LENGTH = "fl_inv_bucket_length";
    private static final String LIMITED_GENE_IDS = "limited_gene_ids"; // for testing
    private static final String LIMITED_CHROMOSOMES = "limited_chromosomes"; // for testing

    private static final String GENE_PAIR_DELIM = "_";

    public static int PRE_GENE_3P_DISTANCE = 10000;
    public static int SHORT_INV_BUCKET = 100000;
    public static long MIN_BUCKET_LENGTH = 100;

    private static final Logger LOGGER = LogManager.getLogger(SvGeneTranscriptCollection.class);

    public FusionLikelihood()
    {
        mDelBucketLengths = Lists.newArrayList();
        mDupBucketLengths = Lists.newArrayList();
        mChrForwardGeneDataMap = Maps.newHashMap();
        mChrReverseGeneDataMap = Maps.newHashMap();
        mDelGenePairCounts = Maps.newHashMap();
        mDupGenePairCounts = Maps.newHashMap();
        mGeneRangeDataMap = Maps.newHashMap();
        mChrPhaseRegionsPosStrand = Maps.newHashMap();
        mChrPhaseRegionsNegStrand = Maps.newHashMap();
        mRestrictedChromosomes = Lists.newArrayList();
        mRestrictedGeneIds = Lists.newArrayList();
        mArmLengthFactor = 0;
        mLogVerbose = false;
    }

    public final Map<String, List<GeneRangeData>> getChrGeneRangeDataMap() { return mChrForwardGeneDataMap; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(DEL_BUCKET_LENGTHS, true, "Semi-colon separated DEL bucket lengths");
        options.addOption(DUP_BUCKET_LENGTHS, true, "Semi-colon separated DUP bucket lengths");
        options.addOption(SHORT_INV_BUCKET_LENGTH, true, "INV bucket length");
        options.addOption(LIMITED_GENE_IDS, true, "List of geneIds to test with");
        options.addOption(LIMITED_CHROMOSOMES, true, "List of chromosomes to test with");
    }

    public void initialise(final CommandLine cmdLineArgs, final SvGeneTranscriptCollection geneTransCache)
    {
        mGeneTransCache = geneTransCache;

        if(cmdLineArgs.hasOption(DEL_BUCKET_LENGTHS))
        {
            setBucketLengths(cmdLineArgs.getOptionValue(DEL_BUCKET_LENGTHS), mDelBucketLengths);
        }

        if(cmdLineArgs.hasOption(DUP_BUCKET_LENGTHS))
        {
            setBucketLengths(cmdLineArgs.getOptionValue(DUP_BUCKET_LENGTHS), mDupBucketLengths);
        }

        if(cmdLineArgs.hasOption(SHORT_INV_BUCKET_LENGTH))
        {
            SHORT_INV_BUCKET = Integer.parseInt(cmdLineArgs.getOptionValue(SHORT_INV_BUCKET_LENGTH));
        }

        if(cmdLineArgs.hasOption(LIMITED_CHROMOSOMES))
        {
            mRestrictedChromosomes = Arrays.stream(cmdLineArgs.getOptionValue(LIMITED_CHROMOSOMES).split(";"))
                    .collect(Collectors.toList());
        }

        // sum up all arm lengths to adjusted same-arm fusion rates
        long maxBucketLength = max(getMaxBucketLength(), 4000000);

        for(Map.Entry<String,Integer> entry : CHROMOSOME_LENGTHS.entrySet())
        {
            final String chromosome = entry.getKey();

            long armLength = getChromosomalArmLength(chromosome, CHROMOSOME_ARM_P);

            if(armLength > maxBucketLength)
                mArmLengthFactor += pow(armLength - maxBucketLength, 2);

            armLength = getChromosomalArmLength(chromosome, CHROMOSOME_ARM_Q);

            if(armLength > maxBucketLength)
                mArmLengthFactor += pow(armLength - maxBucketLength, 2);
        }
    }

    public void setRestrictedGeneIds(final List<String> geneIds) { mRestrictedGeneIds.addAll(geneIds); }

    @VisibleForTesting
    public void initialise(final SvGeneTranscriptCollection geneTransCache, final List<Long> delLengths, final List<Long> dupLengths,
            int shortInvBucketLength, int preGeneDistance)
    {
        mGeneTransCache = geneTransCache;
        mDelBucketLengths.addAll(delLengths);
        mDupBucketLengths.addAll(dupLengths);
        SHORT_INV_BUCKET = shortInvBucketLength;
        PRE_GENE_3P_DISTANCE = preGeneDistance;
    }

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }

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

    public void generateGenePhasingCounts()
    {
        // for each gene, walk through all its transcripts and count up the number of bases in each phasing region (eg 5'UTR, 0-2)
        // also convert these into combined phasing regions (eg where say phase 0 and 1 overlap) for using in downstream analysis
        // finally gather up the total bases by phasing type across each chromosome
        final Map<String, List<EnsemblGeneData>> chrGeneDataMap = mGeneTransCache.getChrGeneDataMap();
        final Map<String, List<EnsemblGeneData>> chrReverseGeneDataMap = mGeneTransCache.getChrReverseGeneDataMap();

        for(Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            if(!mRestrictedChromosomes.isEmpty() && !mRestrictedChromosomes.contains(chromosome))
                continue;

            LOGGER.info("generating phase counts for chromosome({})", chromosome);

            final List<EnsemblGeneData> geneDataList = entry.getValue();
            final List<EnsemblGeneData> reverseGeneDataList = chrReverseGeneDataMap.get(chromosome);

            List<GeneRangeData> geneList = Lists.newArrayList();
            List<GeneRangeData> geneEndFirstList = Lists.newArrayList();
            List<GenePhaseRegion> chrPhaseRegionsPosStrand = Lists.newArrayList();
            List<GenePhaseRegion> chrPhaseRegionsNegStrand = Lists.newArrayList();
            long posStrandTotal = 0;
            long negStrandTotal = 0;

            for(final EnsemblGeneData geneData :entry.getValue())
            {
                GeneRangeData geneRangeData = new GeneRangeData(geneData);

                geneList.add(geneRangeData);

                int index = 0;
                for(; index < geneEndFirstList.size(); ++index)
                {
                    final GeneRangeData rgd = geneEndFirstList.get(index);

                    if(geneData.GeneEnd < rgd.GeneData.GeneEnd)
                        break;
                }

                geneEndFirstList.add(index, geneRangeData);

                if(!mRestrictedGeneIds.isEmpty() && !mRestrictedGeneIds.contains(geneData.GeneId))
                    continue;

                List<GenePhaseRegion> nonOverlappingRegions = null;

                if(!mGeneRangeDataMap.isEmpty())
                {
                    final String phaseData = mGeneRangeDataMap.get(geneData.GeneId);

                    if(phaseData == null || phaseData.isEmpty())
                        continue;

                    geneRangeData.loadRegionsFromCsv(phaseData);
                    nonOverlappingRegions = geneRangeData.getCombinedPhaseRegions();
                }

                if(nonOverlappingRegions == null)
                {
                    // load from Ensembl transcript and exon data
                    final List<TranscriptExonData> transExonDataList = mGeneTransCache.getTransExonData(geneData.GeneId);

                    if (transExonDataList == null)
                        continue;

                    long precedingGeneSAPos = mGeneTransCache.findPrecedingGeneSpliceAcceptorPosition(
                            geneData, geneData.Strand == 1 ? reverseGeneDataList : geneDataList);

                    List<GenePhaseRegion> phaseRegions = generateGenePhaseRegions(geneData, transExonDataList, precedingGeneSAPos);
                    geneRangeData.addPhaseRegions(phaseRegions);

                    // convert to a non-overlapping combined set of phase regions
                    nonOverlappingRegions = Lists.newArrayList();

                    for (GenePhaseRegion region : geneRangeData.getPhaseRegions())
                    {
                        if (region.length() <= 0)
                        {
                            LOGGER.warn("gene({}) invalid length({})", region.GeneId, region.length());
                            continue;
                        }

                        checkAddCombinedGenePhaseRegion(region, nonOverlappingRegions);
                    }

                    geneRangeData.setCombinedPhaseRegions(nonOverlappingRegions);
                }

                if(mLogVerbose)
                {
                    for(GenePhaseRegion region : geneRangeData.getPhaseRegions())
                    {
                        LOGGER.debug("gene({}) range({} -> {}) phase({}) preGene({})",
                                region.GeneId, region.start(), region.end(), region.Phase, region.isAnyPreGene());
                    }
                }

                for(GenePhaseRegion region : nonOverlappingRegions)
                {
                    // validity check
                    if(region.length() <= 0)
                    {
                        LOGGER.error("invalid region: gene({}) range({} -> {}) phase({})",
                                region.GeneId, region.start(), region.end(), region.getCombinedPhase());
                        continue;
                    }

                    // add to the chromosome's summary phasing regions by strand, without concern for overlaps
                    // across genes on the same strand
                    if(geneData.Strand == 1)
                    {
                        posStrandTotal += region.length();
                        addPhaseRegion(chrPhaseRegionsPosStrand, region);
                    }
                    else
                    {
                        negStrandTotal += region.length();
                        addPhaseRegion(chrPhaseRegionsNegStrand, region);
                    }
                }
            }

            mChrPhaseRegionsPosStrand.put(chromosome, chrPhaseRegionsPosStrand);
            mChrPhaseRegionsNegStrand.put(chromosome, chrPhaseRegionsNegStrand);

            LOGGER.debug("chr({}) posStrandTotal({}) negStrandTotal({})", chromosome, posStrandTotal, negStrandTotal);

            mChrForwardGeneDataMap.put(chromosome, geneList);
            mChrReverseGeneDataMap.put(chromosome, geneEndFirstList);
        }
    }

    private void addPhaseRegion(List<GenePhaseRegion> regions, final GenePhaseRegion newRegion)
    {
        // add if new, else expand existing
        for(GenePhaseRegion region : regions)
        {
            if(region.getCombinedPhase() == newRegion.getCombinedPhase())
            {
                region.setEnd(region.end() + newRegion.length());
                return;
            }
        }

        regions.add(new GenePhaseRegion("", 0, newRegion.length(),
                newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus()));
    }

    public static List<GenePhaseRegion> generateGenePhaseRegions(
            final EnsemblGeneData geneData, final List<TranscriptExonData> transExonDataList, long precedingGeneSAPos)
    {
        // 5-prime rules - must be post-promotor
        // 3-prime rules: must be coding and > 1 exon, needs to find first splice acceptor and then uses its phasing
        // need to mark coding and non-coding regions
        List<GenePhaseRegion> phaseRegions = Lists.newArrayList();

        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

        while (!transcriptExons.isEmpty())
        {
            for (int i = 0; i < transcriptExons.size() - 1; ++i)
            {
                // skip single-exon transcripts since without an intronic section their fusion likelihood is negligible
                if(transcriptExons.size() == 1)
                    break;

                TranscriptExonData exonData = transcriptExons.get(i);

                if (exonData.CodingStart == null)
                {
                    // mark the whole transcript as a single UTR
                    GenePhaseRegion phaseRegion = new GenePhaseRegion(
                            geneData.GeneId, exonData.TransStart, exonData.TransEnd, PHASE_NON_CODING);

                    checkAddGenePhaseRegion(phaseRegion, phaseRegions);
                    break;
                }

                TranscriptExonData nextExonData = transcriptExons.get(i+1);

                if(geneData.Strand == 1 && nextExonData.ExonStart > exonData.CodingEnd)
                    break;
                else if(geneData.Strand == -1 && exonData.ExonEnd < exonData.CodingStart)
                    continue;

                // add an upstream gene region (only applicable for downstream fusion genes)
                // with a phase of -1 unless coding starts in the first exon or on the first base of the second exon
                if(precedingGeneSAPos > 0 && geneData.Strand == 1 && i == 0)
                {
                    int regionPhase = (exonData.CodingStart < exonData.ExonEnd) ? nextExonData.ExonPhase : -1;

                    long preDistance = max(exonData.TransStart - precedingGeneSAPos, 0);
                    long upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                    long regionStart = exonData.TransStart - upstreamDistance;
                    long regionEnd = exonData.ExonStart - 1;

                    if(regionStart < regionEnd)
                    {
                        GenePhaseRegion phaseRegion =
                                new GenePhaseRegion(geneData.GeneId, regionStart, regionEnd, mapExonPhase(regionPhase));
                        phaseRegion.setPreGene(true, phaseRegion.Phase);

                        checkAddGenePhaseRegion(phaseRegion, phaseRegions);
                    }
                }
                else if(precedingGeneSAPos > 0 && geneData.Strand == -1 && i == transcriptExons.size() - 2)
                {
                    int regionPhase = (exonData.CodingEnd > nextExonData.ExonStart) ? exonData.ExonPhase : -1;

                    long regionStart = nextExonData.ExonEnd + 1;

                    long preDistance = max(precedingGeneSAPos - exonData.TransEnd, 0);
                    long upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                    long regionEnd = nextExonData.TransEnd + upstreamDistance;

                    if(regionStart < regionEnd)
                    {
                        GenePhaseRegion phaseRegion =
                                new GenePhaseRegion(geneData.GeneId, regionStart, regionEnd, mapExonPhase(regionPhase));
                        phaseRegion.setPreGene(true, phaseRegion.Phase);

                        checkAddGenePhaseRegion(phaseRegion, phaseRegions);
                    }
                }

                // turn the intronic section into a phase region (and fold the exon in with same phasing for simplicity

                GenePhaseRegion phaseRegion = null;

                if(geneData.Strand == 1)
                {
                    phaseRegion = new GenePhaseRegion(
                            geneData.GeneId, exonData.ExonStart, nextExonData.ExonStart, mapExonPhase(exonData.ExonPhaseEnd));
                }
                else
                {
                    phaseRegion = new GenePhaseRegion(
                            geneData.GeneId, exonData.ExonEnd, nextExonData.ExonEnd, mapExonPhase(nextExonData.ExonPhaseEnd));
                }

                checkAddGenePhaseRegion(phaseRegion, phaseRegions);
            }

            teIndex += transcriptExons.size();
            transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
        }

        return phaseRegions;
    }

    public void generateNonProximateCounts()
    {
        // for each arm and each gene in that arm sum up the overlapping counts against all genes beyond the specified
        // DEL and DUP max bucket length, and then all overlapping counts on remote arms
        long proximateLimit = getMaxBucketLength();

        for (Map.Entry<String, List<GeneRangeData>> entry : mChrForwardGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            LOGGER.info("calculating local and remote overlap counts for chr({})", chromosome);

            List<GeneRangeData> geneList = entry.getValue();

            for(int i = 0; i < geneList.size(); ++i)
            {
                GeneRangeData gene1 = geneList.get(i);

                /*
                if(gene1.GeneData.GeneId.equals("ENSG00000189283"))
                {
                    LOGGER.info("specific gene: {}", gene1.GeneData.GeneId);
                }
                */

                // determine the remote overlap count by comparing this gene's regions to all matching remote ones
                for(int j = 0; j <= 1; ++j)
                {
                    final Map<String, List<GenePhaseRegion>> chrPhaseRegions =
                            (j == 0) ? mChrPhaseRegionsPosStrand : mChrPhaseRegionsNegStrand;

                    for (Map.Entry<String, List<GenePhaseRegion>> chrEntry : chrPhaseRegions.entrySet())
                    {
                        if (chrEntry.getKey().equals(chromosome)) // skip the same chromosome
                            continue;

                        final List<GenePhaseRegion> remoteRegions = chrEntry.getValue();

                        for (GenePhaseRegion region : gene1.getCombinedPhaseRegions())
                        {
                            for (GenePhaseRegion remoteRegion : remoteRegions)
                            {
                                // test gene as an downstream vs all remote upstream phasing regions
                                if (hasAnyPhaseMatch(remoteRegion, region, false) || regionsPhaseMatched(remoteRegion, region))
                                {
                                    gene1.addBaseOverlapCountDownstream(NON_PROX_TYPE_REMOTE, region.length() * remoteRegion.length());
                                }

                                // then as an upstream partner
                                if (hasAnyPhaseMatch(region, remoteRegion, false) || regionsPhaseMatched(region, remoteRegion))
                                {
                                    gene1.addBaseOverlapCountUpstream(NON_PROX_TYPE_REMOTE, region.length() * remoteRegion.length());
                                }
                            }
                        }
                    }
                }

                // and now the local ones outside the DEL and DUP proximity lengths, ignoring any strand checking
                // this is considering fusions between DUPs, DUPs or INVs
                for(int j = i+1; j < geneList.size(); ++j)
                {
                    GeneRangeData gene2 = geneList.get(j);

                    if(!gene1.Arm.equals(gene2.Arm))
                        break;

                    int type = -1;

                    if (abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneStart) > proximateLimit
                    && abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneEnd) > proximateLimit
                    && abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneEnd) > proximateLimit
                    && abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneStart) > proximateLimit)
                    {
                        type = NON_PROX_TYPE_LONG_SAME_ARM;
                    }
                    else if(gene1.GeneData.Strand != gene2.GeneData.Strand
                    && abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneStart) < SHORT_INV_BUCKET
                    && abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneEnd) < SHORT_INV_BUCKET
                    && abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneEnd) < SHORT_INV_BUCKET
                    && abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneStart) < SHORT_INV_BUCKET)
                    {
                        type = NON_PROX_TYPE_SHORT_INV;
                    }
                    else
                    {
                        continue;
                    }

                    // calculate phasing overlap areas
                    for (GenePhaseRegion region1 : gene1.getCombinedPhaseRegions())
                    {
                        // the downstream gene of the potential fusion cannot be non-coding
                        for (GenePhaseRegion region2 : gene2.getCombinedPhaseRegions())
                        {
                            long regionOverlap = region1.length() * region2.length();

                            if(hasAnyPhaseMatch(region1, region2, false))
                            {
                                gene1.addBaseOverlapCountUpstream(type, regionOverlap);
                                gene1.addBaseOverlapCountDownstream(type, regionOverlap);

                                gene2.addBaseOverlapCountUpstream(type, regionOverlap);
                                gene2.addBaseOverlapCountDownstream(type, regionOverlap);
                            }
                            else if(regionsPhaseMatched(region1, region2))
                            {
                                gene1.addBaseOverlapCountUpstream(type, regionOverlap);
                                gene2.addBaseOverlapCountDownstream(type, regionOverlap);
                            }
                            else if(regionsPhaseMatched(region2, region1))
                            {
                                gene1.addBaseOverlapCountDownstream(type, regionOverlap);
                                gene2.addBaseOverlapCountUpstream(type, regionOverlap);
                            }
                        }
                    }
                }
            }
        }
    }

    public void generateProximateFusionCounts()
    {
        if(mDelBucketLengths.isEmpty() && mDupBucketLengths.isEmpty())
            return;

        // walk through each chromosome, and test distances between each gene and the next ones for potential
        // fusions from DELs and DUPs within the configured bucket distances

        LOGGER.info("finding proximate DEL fusion candidates");

        for (Map.Entry<String, List<GeneRangeData>> entry : mChrForwardGeneDataMap.entrySet())
        {
            LOGGER.info("proximate forward-strand fusion candidates from chromosome({})", entry.getKey());
            findProximateFusions(entry.getValue(), 1);
        }

        LOGGER.info("finding proximate DUP fusion candidates");

        for (Map.Entry<String, List<GeneRangeData>> entry : mChrReverseGeneDataMap.entrySet())
        {
            LOGGER.info("proximate reverse-strand fusion candidates from chromosome({})", entry.getKey());
            findProximateFusions(entry.getValue(), -1);
        }
    }

    private long getMaxBucketLength()
    {
        long delLimit = !mDelBucketLengths.isEmpty() ? mDelBucketLengths.get(mDelBucketLengths.size() - 1) : MIN_BUCKET_LENGTH;
        long dupLimit = !mDupBucketLengths.isEmpty() ? mDupBucketLengths.get(mDupBucketLengths.size() - 1) : MIN_BUCKET_LENGTH;
        return max(delLimit, dupLimit);
    }

    private void findProximateFusions(List<GeneRangeData> geneRangeList, int strandMatch)
    {
        long rangeLimit = getMaxBucketLength();

        for(int lowerIndex = 0; lowerIndex < geneRangeList.size(); ++lowerIndex)
        {
            GeneRangeData lowerGene = geneRangeList.get(lowerIndex);

            if(lowerGene.GeneData.Strand != strandMatch)
                continue;

            // allow same-gene fusions, so start the next gene at the current one
            for(int upperIndex = lowerIndex; upperIndex < geneRangeList.size(); ++upperIndex)
            {
                GeneRangeData upperGene = geneRangeList.get(upperIndex);

                if(upperGene.GeneData.Strand != strandMatch)
                    continue;

                // exit if the gene is now too far away from the lower gene
                if(upperGene.GeneData.GeneStart - lowerGene.GeneData.GeneEnd > rangeLimit)
                    break;

                findProximateFusions(lowerGene, upperGene);
            }
        }
    }

    private void findProximateFusions(GeneRangeData lowerGene, GeneRangeData upperGene)
    {
        // find all matching phasing regions and account for any duplicated overlapping regions from different phasing matches
        // for all those found, assign them to a length bucket and find the number of bases that fall into the length bucket region

        boolean isForwardStrand = (lowerGene.GeneData.Strand == 1);

        for(int i = 0; i <= 1; ++i)
        {
            boolean isDel = (i == 0);
            boolean lowerGeneIsUpstream = (isDel == isForwardStrand);
            boolean upperGeneIsUpstream = !lowerGeneIsUpstream;

            // first find mutually exclusive matching regions
            for (GenePhaseRegion lowerRegion : lowerGene.getCombinedPhaseRegions())
            {
                // the downstream gene of the potential fusion cannot be non-coding
                if(!lowerGeneIsUpstream && lowerRegion.hasPhaseOnly(PHASE_NON_CODING))
                    continue;

                for (GenePhaseRegion upperRegion : upperGene.getCombinedPhaseRegions())
                {
                    // ignore overlapping regions for now since it's not clear whether a DUP or DEL would be required
                    if (!(lowerRegion.end() < upperRegion.start() || lowerRegion.start() > upperRegion.end()))
                        continue;

                    // if (lowerRegion.end() > upperRegion.start())
                    //    continue;

                    boolean phaseMatched = hasAnyPhaseMatch(lowerRegion, upperRegion, false);

                    if(!phaseMatched)
                    {
                        if(upperGeneIsUpstream && regionsPhaseMatched(upperRegion, lowerRegion))
                            phaseMatched = true;
                        else if(!upperGeneIsUpstream && regionsPhaseMatched(lowerRegion, upperRegion))
                            phaseMatched = true;
                    }

                    if(!phaseMatched)
                        continue;

                    List<Long> bucketLengths = isDel ? mDelBucketLengths : mDupBucketLengths;

                    Map<Integer, Long> bucketOverlapCounts = calcOverlapBucketAreas(
                            bucketLengths, lowerGene, upperGene, lowerRegion, upperRegion, isDel);

                    if(!bucketOverlapCounts.isEmpty())
                    {
                        if (mLogVerbose)
                        {
                            LOGGER.debug("gene({}: {}) and gene({}: {}) matched-region lower({} -> {}) upper({} -> {}) phase({}):",
                                    lowerGene.GeneData.GeneId, lowerGene.GeneData.GeneName,
                                    upperGene.GeneData.GeneId, upperGene.GeneData.GeneName,
                                    lowerRegion.start(), lowerRegion.end(), upperRegion.start(), upperRegion.end(),
                                    lowerRegion.getCombinedPhase());
                        }

                        for (Map.Entry<Integer, Long> entry : bucketOverlapCounts.entrySet())
                        {
                            int bucketIndex = entry.getKey();
                            long overlap = entry.getValue();
                            addGeneFusionData(lowerGene, upperGene, overlap, isDel, bucketIndex);

                            if(mLogVerbose)
                            {
                                long bucketLen = bucketLengths.get(bucketIndex);
                                LOGGER.debug("matched-region bucketLength({}: {}) overlap({})", bucketIndex, bucketLen, overlap);
                            }
                        }
                    }
                }
            }
        }
    }

    public static Map<Integer, Long> calcOverlapBucketAreas(
            final List<Long> bucketLengths, GeneRangeData lowerGene, GeneRangeData upperGene,
            GenePhaseRegion lowerRegion, GenePhaseRegion upperRegion, boolean isDel)
    {
        Map<Integer, Long> bucketOverlapCounts = Maps.newHashMap();

        if(lowerRegion.length() < 0 || upperRegion.length() < 0)
        {
            LOGGER.warn("negative region lengths");
        }

        for (int i = 0; i < bucketLengths.size() - 1; ++i)
        {
            long minBucketLen = bucketLengths.get(i);
            long maxBucketLen = bucketLengths.get(i + 1);

            // first check whether the bucket can link these 2 genes
            if(lowerRegion.start() + minBucketLen >= upperRegion.end() || lowerRegion.end() + maxBucketLen <= upperRegion.start())
            {
                continue;
            }

            /* Example scenarios:
                - gene ranges 100-200 and 300-400
                - if min bucket length is less than 100 (ie the distance between the lower gene end and upper gene start) then no restriction
                - if max bucket length is greater than 300 (ie the distance from lower gene start to upper gene end) then no restriction
                - in this case the base overlap is 100 x 100, so 10K - the number of different ways a unique SV can fuse these 2 genes
                - if the min bucket length were increased say to 225, then the overlap ranges is limited
                - if the max bucket length were a limiting factor at say 275 then 100 + 275 only reaches 375
                - if the max bucket length = 150 then 100 + 150 doesn't reach the start of the upper gene, so the base overlap region must start from 150+
            */

            long baseOverlapArea = 0;

            if(minBucketLen <= upperRegion.start() - lowerRegion.end() && maxBucketLen >= upperRegion.end() - lowerRegion.start())
            {
                // no restriction on the overlap
                baseOverlapArea = lowerRegion.length() * upperRegion.length();
            }
            else
            {
                long lowerStart, lowerEnd, upperStart;

                long upperEnd = min(upperRegion.end(), lowerRegion.end() + maxBucketLen);

                if(lowerRegion.start() + minBucketLen > upperRegion.start())
                {
                    // min bucket length limits bases in the upper region
                    lowerStart = lowerRegion.start();
                    upperStart = lowerRegion.start() + minBucketLen;

                    if(lowerRegion.end() + minBucketLen > upperRegion.end())
                    {
                        lowerEnd = upperRegion.end() - minBucketLen;
                    }
                    else
                    {
                        lowerEnd = lowerRegion.end();
                    }
                }
                else
                {
                    // no min length restrictions, so just check if the max length is restrictive
                    upperStart = upperRegion.start();
                    lowerEnd = lowerRegion.end();

                    if(lowerRegion.start() + maxBucketLen < upperRegion.start())
                    {
                        lowerStart = upperRegion.start() - maxBucketLen;
                    }
                    else
                    {
                        lowerStart = lowerRegion.start();
                    }
                }

                for (long base = lowerStart; base <= lowerEnd; ++base)
                {
                    long overlap = min(upperEnd, base + maxBucketLen) - max(upperStart, base + minBucketLen);

                    if(overlap > 0)
                        baseOverlapArea += overlap;

                    // is there any early exit here once overlap starts to be negative?
                }
            }

            setBucketLengthData(isDel ? lowerGene.getDelFusionBaseCounts() : lowerGene.getDupFusionBaseCounts(), i, baseOverlapArea);

            if(!upperGene.GeneData.GeneId.equals(lowerGene.GeneData.GeneId))
            {
                // avoid double counting same gene fusions
                setBucketLengthData(isDel ? upperGene.getDelFusionBaseCounts() : upperGene.getDupFusionBaseCounts(), i, baseOverlapArea);
            }

            bucketOverlapCounts.put(i, baseOverlapArea);
        }

        return bucketOverlapCounts;

    }

    private static void setBucketLengthData(Map<Integer,Long> countsData, int bucketIndex, long newCounts)
    {
        // initialise the array if empty
        Long bucketCount = countsData.get(bucketIndex);

        if(bucketCount == null)
        {
            countsData.put(bucketIndex, newCounts);
        }
        else
        {
            countsData.put(bucketIndex, bucketCount + newCounts);
        }
    }

    private static final int BUCKET_MIN = 0;
    private static final int BUCKET_MAX = 1;

    private long[] getBucketLengthMinMax(boolean isDel, int bucketIndex)
    {
        List<Long> bucketLengths = isDel ? mDelBucketLengths : mDupBucketLengths;

        if(bucketIndex >= bucketLengths.size())
            return new long[2];

        return new long[] {bucketLengths.get(bucketIndex), bucketLengths.get(bucketIndex + 1)};
    }

    public void addGeneFusionData(final GeneRangeData lowerGene, final GeneRangeData upperGene, long overlapCount, boolean isDel, int bucketIndex)
    {
        long[] bucketMinMax = getBucketLengthMinMax(isDel, bucketIndex);
        long bucketWidth = bucketMinMax[BUCKET_MAX] - bucketMinMax[BUCKET_MIN];
        double fusionRate = overlapCount / (bucketWidth * GENOME_BASE_COUNT);

        if(fusionRate < MIN_FUSION_RATE)
            return;

        if(mLogVerbose)
        {
            LOGGER.debug("gene pair({} & {}) adding {} overlap({}) for bucket length index({})",
                    lowerGene.GeneData.GeneName, upperGene.GeneData.GeneName,
                    isDel ? "DEL" : "DUP", overlapCount, bucketIndex);
        }

        Map<String, Map<Integer,Long>> genePairCounts = isDel ? mDelGenePairCounts : mDupGenePairCounts;

        final String genePair = lowerGene.GeneData.GeneId + GENE_PAIR_DELIM + upperGene.GeneData.GeneId;
        Map<Integer,Long> bucketOverlapCounts = genePairCounts.get(genePair);

        if(bucketOverlapCounts == null)
        {
            bucketOverlapCounts = Maps.newHashMap();
            genePairCounts.put(genePair, bucketOverlapCounts);
        }

        setBucketLengthData(bucketOverlapCounts, bucketIndex, overlapCount);
    }

    public static void checkAddGenePhaseRegion(final GenePhaseRegion newRegion, final List<GenePhaseRegion> regions)
    {
        if(newRegion.length() <= 0)
        {
            LOGGER.warn("gene({}) attempting to add invalid region with length({})", newRegion.GeneId, newRegion.length());
            return;
        }

        // add a region if it doesn't overlap with any existing
        // otherwise expand the existing region
        for(final GenePhaseRegion region : regions)
        {
            if(region.Phase != newRegion.Phase)
                continue;

            if(region.start() > newRegion.end() || region.end() < newRegion.start())
                continue;

            // widen the region to cover both region boundaries
            region.setStart(min(region.start(), newRegion.start()));
            region.setEnd(max(region.end(), newRegion.end()));

            if(region.length() <= 0)
            {
                LOGGER.warn("gene({}) adjusted to invalid region with length({})", region.GeneId, region.length());
                return;
            }

            return;
        }

        // add new non-overlapping or differently-phased region
        regions.add(newRegion);
    }

    public static void checkAddCombinedGenePhaseRegion(final GenePhaseRegion regionToAdd, final List<GenePhaseRegion> regions)
    {
        // combine all regions into non-overlapping regions, allowing mixed phasings to do so

        // clone to avoid changing the contents of the original
        GenePhaseRegion tmp = new GenePhaseRegion(regionToAdd.GeneId, regionToAdd.start(), regionToAdd.end(),
                regionToAdd.getPhaseArray(), regionToAdd.getPreGenePhaseStatus());

        List<GenePhaseRegion> newRegions = Lists.newArrayList(tmp);

        // look for overlapping regions and combine or split them as required
        // split any new regions until they can be added without any further splits
        while(!newRegions.isEmpty())
        {
            GenePhaseRegion newRegion = newRegions.get(0);
            newRegions.remove(0);

            int index = 0;
            boolean regionSplit = false;
            while(index < regions.size())
            {
                GenePhaseRegion region = regions.get(index);

                if(newRegion == null || region == null)
                {
                    LOGGER.error("isNull: region({}) or new region({})", region == null, newRegion == null);
                    return;
                }

                if (region.start() > newRegion.end() || region.end() < newRegion.start())
                {
                    ++index;
                    continue;
                }

                // preserve pre-gene region status only if both regions being combined are pre-gene, otherwise the region will
                // function like a normal phasing region

                // if the region matches on combined phase exactly, then expand the new region to cover both, remove the existing and continue
                if (region.getCombinedPhase() == newRegion.getCombinedPhase())
                {
                    newRegion.setStart(min(region.start(), newRegion.start()));
                    newRegion.setEnd(max(region.end(), newRegion.end()));
                    newRegion.addPhases(region.getPhaseArray(), region.getPreGenePhaseStatus());
                    regions.remove(index);
                    continue;
                }

                // otherwise split and combine the regions so there are no overlaps
                regionSplit = true;

                // first check for one region enclosing another, the overlaps
                GenePhaseRegion extraRegion = null;

                if (region.start() == newRegion.start() && region.end() == newRegion.end())
                {
                    region.addPhases(newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());
                    newRegion.setEnd(newRegion.start()); // won't be added
                }
                else if(region.end() == newRegion.start() || region.start() == newRegion.end())
                {
                    // a precise overlap move base 1 apart and continue
                    if(region.end() == newRegion.start())
                    {
                        region.setEnd(newRegion.start() - 1);
                    }
                    else
                    {
                        newRegion.setEnd(region.start() - 1);
                    }
                }
                else if (region.start() <= newRegion.start() && region.end() >= newRegion.end())
                {
                    // existing region enclosed the new region
                    // split the outer region in 2 and make a combined inner region
                    extraRegion = new GenePhaseRegion(region.GeneId, newRegion.end() + 1, region.end(),
                            region.getPhaseArray(), region.getPreGenePhaseStatus());

                    newRegion.addPhases(region.getPhaseArray(), region.getPreGenePhaseStatus());

                    region.setEnd(newRegion.start() - 1);
                }
                else if (newRegion.start() <= region.start() && newRegion.end() >= region.end())
                {
                    // existing region falls within new region
                    extraRegion = new GenePhaseRegion(newRegion.GeneId, region.end() + 1, newRegion.end(),
                            newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());

                    region.addPhases(newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());

                    newRegion.setEnd(region.start() - 1);

                }
                else if (newRegion.start() <= region.start())
                {
                    // new region precedes and overlaps the existing

                    // create a new region for the overlapping section
                    extraRegion = new GenePhaseRegion(newRegion.GeneId, region.start(), newRegion.end(),
                            newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());

                    extraRegion.addPhases(region.getPhaseArray(), region.getPreGenePhaseStatus());

                    long regionStart = newRegion.end() + 1;
                    long newRegionEnd = region.start() - 1;

                    region.setStart(regionStart);

                    newRegion.setEnd(newRegionEnd);
                }
                else if (region.start() <= newRegion.start())
                {
                    // existing region precedes the new region
                    extraRegion = new GenePhaseRegion(newRegion.GeneId, newRegion.start(), region.end(),
                            newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());

                    extraRegion.addPhases(region.getPhaseArray(), region.getPreGenePhaseStatus());

                    long regionEnd = newRegion.start() - 1;
                    long newRegionStart = region.end() + 1;

                    region.setEnd(regionEnd);

                    newRegion.setStart(newRegionStart);
                }

                if(newRegion != null && newRegion.length() >= 1)
                {
                    newRegions.add(newRegion);
                }

                if(extraRegion != null && extraRegion.length() >= 1)
                {
                    newRegions.add(extraRegion);
                }

                ++index;
            }

            if(!regionSplit)
            {
                regions.add(newRegion);
            }
        }

        // remove 0-1 bases regions
        int index = 0;
        while(index < regions.size())
        {
            if(regions.get(index).length() <= 1)
            {
                regions.remove(index);
            }
            else
            {
                ++index;
            }
        }
    }

    public void writeGeneLikelihoodData(final String outputDir)
    {
        LOGGER.info("writing output files");

        writeGeneData(outputDir);
        writeTopProximateFusionCandidates(outputDir);

        if(mGeneRangeDataMap.isEmpty())
            writeGenePhaseData(outputDir);
    }

    private static final String GENE_PHASE_DATA_FILE = "GFL_GENE_PHASE_DATA.csv";

    private void writeGenePhaseData(final String outputDir)
    {
        LOGGER.info("writing gene phase data cache file");

        try
        {
            String outputFilename = outputDir + GENE_PHASE_DATA_FILE;

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("GeneId,PhaseRegions");
            writer.newLine();

            for (Map.Entry<String, List<GeneRangeData>> entry : mChrForwardGeneDataMap.entrySet())
            {
                for (final GeneRangeData geneData : entry.getValue())
                {
                    if(geneData.getPhaseRegions().isEmpty())
                        continue;

                    writer.write(geneData.toCsv());
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

    public boolean loadGenePhaseData(final String inputDir)
    {
        // attempt to load cached gene phase data from file to avoid re-creating it each time
        String filename = inputDir + GENE_PHASE_DATA_FILE;

        if(!Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                return false;
            }

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(PGD_DELIMITER);

                if(items.length != 3)
                    continue;

                final String geneId = items[0];
                mGeneRangeDataMap.put(geneId, items[1] + PGD_DELIMITER + items[2]);
            }
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read gene range data CSV file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    private void writeGeneData(final String outputDir)
    {
        try
        {
            String outputFilename = outputDir + "GFL_GENE_DATA.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("GeneId,GeneName,Chromosome,Arm,GeneStart,GeneEnd,Strand");
            writer.write(",FivePrimeUTR,Phase0,Phase1,Phase2,NonCoding,PreGene");
            writer.write(",ShortInvRateUp,ShortInvRateDown,SameArmRateUp,SameArmRateDown,RemoteRateUp,RemoteRateDown");
            writer.newLine();

            // adjustment factors to convert overlap base count into rates
            double remoteFusionFactor = 1.0 / (GENOME_BASE_COUNT * GENOME_BASE_COUNT);
            double sameArmFusionFactor = 1.0 / mArmLengthFactor;
            double shortInvFusionFactor = 1.0 / (SHORT_INV_BUCKET * GENOME_BASE_COUNT);

            for(Map.Entry<String, List<GeneRangeData>> entry : mChrForwardGeneDataMap.entrySet())
            {
                for(final GeneRangeData geneData :entry.getValue())
                {
                    long phase0 = 0;
                    long phase1 = 0;
                    long phase2 = 0;
                    long phase5pUTR = 0;
                    long nonCoding = 0;
                    long preGene = 0;

                    for(final GenePhaseRegion region : geneData.getPhaseRegions())
                    {
                        if(region.Phase == PHASE_5P_UTR)
                            phase5pUTR += region.length();
                        else if(region.Phase == PHASE_NON_CODING)
                            nonCoding += region.length();
                        else if(region.Phase == PHASE_0)
                            phase0 += region.length();
                        else if(region.Phase == PHASE_1)
                            phase1 += region.length();
                        else if(region.Phase == PHASE_2)
                            phase2 += region.length();

                        if(region.isAnyPreGene())
                            preGene += region.length();
                    }

                    writer.write(String.format("%s,%s,%s,%s,%d,%d,%d",
                            geneData.GeneData.GeneId, geneData.GeneData.GeneName, geneData.GeneData.Chromosome, geneData.Arm,
                            geneData.GeneData.GeneStart, geneData.GeneData.GeneEnd, geneData.GeneData.Strand));

                    writer.write(String.format(",%d,%d,%d,%d,%d,%d",
                            phase5pUTR, phase0, phase1, phase2, nonCoding, preGene));

                    writer.write(String.format(",%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                            geneData.getBaseOverlapCountUpstream(NON_PROX_TYPE_SHORT_INV) * shortInvFusionFactor,
                            geneData.getBaseOverlapCountDownstream(NON_PROX_TYPE_SHORT_INV) * shortInvFusionFactor,
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

    private static final double GENOME_BASE_COUNT = 3e9;
    private static final double MIN_FUSION_RATE = 1e-12;

    private void writeTopProximateFusionCandidates(final String outputDir)
    {
        LOGGER.info("total gene-pair candidate count: dels({}) dups({})", mDelGenePairCounts.size(), mDupGenePairCounts.size());

        try
        {
            String outputFilename = outputDir + "GFL_DEL_DUP_PROXIMATES.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("Type,LengthMin,LengthMax,GeneIdUp,GeneNameUp,GeneIdDown,GeneNameDown,Chromosome,Strand,ProximateRate");
            writer.newLine();

            for(int i = 0; i <= 1; ++i)
            {
                boolean isDel = (i == 0);
                Map<String, Map<Integer, Long>> genePairCounts = isDel ? mDelGenePairCounts : mDupGenePairCounts;

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

                        long[] bucketMinMax = getBucketLengthMinMax(isDel, bucketIndex);
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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addCmdLineArgs(options);
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log in verbose mode");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Ensembl gene transcript data cache directory");

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

        boolean hasCachedData = fusionLikelihood.loadGenePhaseData(outputDir);

        boolean limitedLoading = !restrictedGeneIds.isEmpty();

        if(!ensemblDataCache.loadEnsemblData(limitedLoading || hasCachedData))
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

        if(restrictedGeneIds.size() == 1)
            fusionLikelihood.setLogVerbose(true);

        fusionLikelihood.generateGenePhasingCounts();
        fusionLikelihood.generateProximateFusionCounts();
        fusionLikelihood.generateNonProximateCounts();
        fusionLikelihood.writeGeneLikelihoodData(outputDir);

        LOGGER.info("Gene likelihood data generation complete");
    }
}
