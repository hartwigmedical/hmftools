package com.hartwig.hmftools.svanalysis.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svanalysis.fusion.GeneRangeData.GENE_PHASING_REGION_5P_UTR;
import static com.hartwig.hmftools.svanalysis.fusion.GeneRangeData.GENE_PHASING_REGION_CODING_0;
import static com.hartwig.hmftools.svanalysis.fusion.GeneRangeData.GENE_PHASING_REGION_CODING_1;
import static com.hartwig.hmftools.svanalysis.fusion.GeneRangeData.GENE_PHASING_REGION_CODING_2;
import static com.hartwig.hmftools.svanalysis.fusion.GeneRangeData.GENE_PHASING_REGION_MAX;
import static com.hartwig.hmftools.svanalysis.fusion.GeneRangeData.GENE_PHASING_REGION_PROMOTOR;
import static com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection.nextTranscriptExons;
import static com.hartwig.hmftools.svanalysis.types.SvaConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.svanalysis.types.SvaConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.svanalysis.types.SvaConfig.LOG_DEBUG;
import static com.hartwig.hmftools.svanalysis.types.SvaConfig.formOutputPath;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.svanalysis.analysis.SvUtilities;
import com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection;

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

    private List<Long> mDelBucketLengths;
    private List<Long> mDupBucketLengths;

    private final Map<String, List<GeneRangeData>> mChrForwardGeneDataMap;
    private final Map<String, List<GeneRangeData>> mChrReverseGeneDataMap;

    private static final String DEL_BUCKET_LENGTHS = "fl_del_bucket_lengths";
    private static final String DUP_BUCKET_LENGTHS = "fl_dup_bucket_lengths";
    private static final String LIMITED_GENE_IDS = "limited_gene_ids"; // for testing

    private static final Logger LOGGER = LogManager.getLogger(SvGeneTranscriptCollection.class);

    public FusionLikelihood()
    {
        mDelBucketLengths = Lists.newArrayList();
        mDupBucketLengths = Lists.newArrayList();
        mChrForwardGeneDataMap = Maps.newHashMap();
        mChrReverseGeneDataMap = Maps.newHashMap();
    }

    public final Map<String, List<GeneRangeData>> getChrGeneRangeDataMap() { return mChrForwardGeneDataMap; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(DEL_BUCKET_LENGTHS, true, "Semi-colon separated DEL bucket lengths");
        options.addOption(DUP_BUCKET_LENGTHS, true, "Semi-colon separated DUP bucket lengths");
        options.addOption(LIMITED_GENE_IDS, true, "List of geneIds to test with");
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
        for(Map.Entry<String, List<EnsemblGeneData>> entry : mGeneTransCache.getChrGeneDataMap().entrySet())
        {
            final String chromosome = entry.getKey();

            List<GeneRangeData> geneList = Lists.newArrayList();
            List<GeneRangeData> geneEndFirstList = Lists.newArrayList();

            for(final EnsemblGeneData geneData :entry.getValue())
            {
                final List<TranscriptExonData> transExonDataList = mGeneTransCache.getTransExonData(geneData.GeneId);

                if(transExonDataList == null)
                    continue;

                GeneRangeData geneRangeData = new GeneRangeData(geneData);

                List<GenePhaseRegion> phaseRegions = generateGenePhaseRegions(geneData, transExonDataList);
                geneRangeData.addPhaseRegions(phaseRegions);

                /*

                setGenePhasingCounts(
                        geneData, transExonDataList,
                        geneRangeData.getPhaseCounts(true),
                        geneRangeData.getPhaseCounts(false));

                */

                geneList.add(geneRangeData);

                int index = 0;
                for(; index < geneEndFirstList.size(); ++index)
                {
                    final GeneRangeData rgd = geneEndFirstList.get(index);

                    if(geneData.GeneEnd < rgd.GeneData.GeneEnd)
                        break;
                }

                geneEndFirstList.add(index, geneRangeData);
            }

            mChrForwardGeneDataMap.put(chromosome, geneList);
            mChrReverseGeneDataMap.put(chromosome, geneEndFirstList);
        }
    }

    private void checkAddGenePhaseRegion(final GenePhaseRegion newRegion, final List<GenePhaseRegion> regions)
    {
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
            return;
        }

        // add new non-overlapping or differently-phased region
        regions.add(newRegion);
    }

    private List<GenePhaseRegion> generateGenePhaseRegions(final EnsemblGeneData geneData, final List<TranscriptExonData> transExonDataList)
    {
        List<GenePhaseRegion> phaseRegions = Lists.newArrayList();

        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

        while (!transcriptExons.isEmpty())
        {
            for (int i = 0; i < transcriptExons.size() - 1; ++i)
            {
                TranscriptExonData exonData = transcriptExons.get(i);

                if (exonData.CodingStart == null)
                {
                    // mark the whole transcript as a single UTR
                    GenePhaseRegion phaseRegion = new GenePhaseRegion(
                            geneData.GeneId, exonData.TransStart, exonData.TransEnd, -1);

                    checkAddGenePhaseRegion(phaseRegion, phaseRegions);
                    break;
                }

                if(geneData.Strand == 1 && exonData.ExonStart > exonData.CodingEnd)
                    break;
                else if(geneData.Strand == -1 && exonData.ExonEnd < exonData.CodingStart)
                    continue;

                TranscriptExonData nextExonData = transcriptExons.get(i+1);

                // turn the intronic section into a phase region (and fold the exon in with same phasing for simplicity

                GenePhaseRegion phaseRegion = null;

                if(geneData.Strand == 1)
                {
                    phaseRegion = new GenePhaseRegion(geneData.GeneId, exonData.ExonStart, nextExonData.ExonStart, exonData.ExonPhaseEnd);
                }
                else
                {
                    phaseRegion = new GenePhaseRegion(geneData.GeneId, exonData.ExonEnd, nextExonData.ExonEnd, nextExonData.ExonPhaseEnd);
                }

                checkAddGenePhaseRegion(phaseRegion, phaseRegions);
            }

            teIndex += transcriptExons.size();
            transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
        }

        return phaseRegions;
    }

    @NotNull
    public static void setGenePhasingCounts(
            final EnsemblGeneData geneData, final List<TranscriptExonData> transExonDataList,
            int[] fivePrimePhaseCounts, int[] threePrimePhaseCounts)
    {
        // sum up the number of bases in each phasing region across all transcripts for the gene
        // split by considering the gene independently as a 5 or 3 prime partner in a fusion

        // 5-prime rules - must be post-promotor (exon 2 onwards)
        // 3-prime rules: must be coding and > 1 exon
        long geneStart = geneData.GeneStart;
        long geneEnd = geneData.GeneEnd;
        int geneLength = (int) (geneEnd - geneStart + 1);

        boolean[][] geneBases5P = new boolean[geneLength][GENE_PHASING_REGION_MAX];
        boolean[][] geneBases3P = new boolean[geneLength][GENE_PHASING_REGION_MAX];

        boolean hasCodingExons = false;
        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

        while (!transcriptExons.isEmpty())
        {
            if (geneData.Strand == 1)
            {
                for (int i = 0; i < transcriptExons.size(); ++i)
                {
                    final TranscriptExonData exonData = transcriptExons.get(i);

                    // non-coding is valid for 5P, not 3P
                    if (exonData.CodingStart == null)
                    {
                        for (long j = exonData.TransStart; j <= exonData.TransEnd; ++j)
                        {
                            int gbPos = (int) (j - geneStart);
                            geneBases5P[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                        }

                        break;
                    }

                    hasCodingExons = true;

                    if (exonData.ExonStart > exonData.CodingEnd) // past end of coding region
                        break;

                    // first mark the entire pre-coding region as 5' UTR
                    if (i == 0)
                    {
                        for (long j = geneStart; j <= exonData.CodingStart; ++j)
                        {
                            int gbPos = (int) (j - geneStart);
                            geneBases5P[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                            geneBases3P[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                        }
                    }

                    if (exonData.ExonEnd < exonData.CodingStart) // already accounted for
                        continue;

                    // now handle the exon's phasing
                    long codingStart = max(exonData.ExonStart, exonData.CodingStart);
                    for (long j = codingStart; j <= exonData.ExonEnd; ++j)
                    {
                        int gbPos = (int) (j - geneStart);

                        if (j > exonData.CodingEnd)
                            break;

                        long adjustedPhase = exonData.ExonPhase + (j - codingStart);
                        int calcPhase = (int) (adjustedPhase % 3);
                        geneBases5P[gbPos][phaseToRegion(calcPhase)] = true;
                        geneBases3P[gbPos][phaseToRegion(calcPhase)] = true;
                    }

                    // fill in the intronic phasing between coding exons
                    if (i < transcriptExons.size() - 1)
                    {
                        final TranscriptExonData nextExon = transcriptExons.get(i + 1);

                        if (nextExon.ExonStart <= nextExon.CodingEnd)
                        {
                            int regionType = phaseToRegion(exonData.ExonPhaseEnd);

                            for (long j = exonData.ExonEnd + 1; j < nextExon.ExonStart; ++j)
                            {
                                int gbPos = (int) (j - geneStart);
                                geneBases5P[gbPos][regionType] = true;
                                geneBases3P[gbPos][regionType] = true;
                            }
                        }
                    }
                }
            }
            else
            {
                // navigate through as per the exon rank
                for (int i = transcriptExons.size() - 1; i >= 0; --i)
                {
                    final TranscriptExonData exonData = transcriptExons.get(i);

                    if (exonData.CodingStart == null)
                    {
                        for (long j = exonData.TransStart; j <= exonData.TransEnd; ++j)
                        {
                            int gbPos = (int) (geneEnd - j);
                            geneBases5P[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                        }

                        break;
                    }

                    hasCodingExons = true;

                    if (exonData.ExonEnd < exonData.CodingStart) // past end of coding region
                        break;

                    // first mark the entire pre-coding region as 5' UTR
                    if (i == transcriptExons.size() - 1)
                    {
                        for (long j = exonData.CodingEnd; j <= geneEnd; ++j)
                        {
                            int gbPos = (int) (geneEnd - j);
                            geneBases5P[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                            geneBases3P[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                        }
                    }

                    if (exonData.ExonStart > exonData.CodingEnd) // already accounted for
                        continue;

                    // allocate the exon's phasing - working backwwards this time
                    long codingEnd = min(exonData.ExonEnd, exonData.CodingEnd);
                    for (long j = codingEnd; j >= exonData.ExonStart; --j)
                    {
                        int gbPos = (int) (geneEnd - j);

                        if (j < exonData.CodingStart)
                            break;

                        long adjustedPhase = exonData.ExonPhase + (codingEnd - j);
                        int calcPhase = (int) (adjustedPhase % 3);
                        geneBases5P[gbPos][phaseToRegion(calcPhase)] = true;
                        geneBases3P[gbPos][phaseToRegion(calcPhase)] = true;
                    }

                    // fill in the intronic phasing between coding exons
                    if (i > 0)
                    {
                        final TranscriptExonData nextExon = transcriptExons.get(i - 1);

                        if (nextExon.ExonEnd >= nextExon.CodingStart)
                        {
                            int regionType = phaseToRegion(exonData.ExonPhaseEnd);

                            for (long j = nextExon.ExonEnd + 1; j < exonData.ExonStart; ++j)
                            {
                                int gbPos = (int) (geneEnd - j);
                                geneBases5P[gbPos][regionType] = true;
                                geneBases3P[gbPos][regionType] = true;
                            }
                        }
                    }
                }
            }

            teIndex += transcriptExons.size();
            transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
        }

        // now compute the number of bases for each phasing region
        for(int i = 0; i < geneLength; ++i)
        {
            for(int j = 0; j < GENE_PHASING_REGION_MAX; ++j)
            {
                if(geneBases5P[i][j])
                    ++fivePrimePhaseCounts[j];

                if(geneBases3P[i][j])
                    ++threePrimePhaseCounts[j];
            }
        }

        if(hasCodingExons)
        {
            LOGGER.debug("gene({}: {}) length({}) region counts: pre-coding({}) 5P phases(0={} 1={} 2={}) 3P phases(0={} 1={} 2={})",
                    geneData.GeneId, geneData.GeneName, geneLength,
                    fivePrimePhaseCounts[GENE_PHASING_REGION_5P_UTR], fivePrimePhaseCounts[GENE_PHASING_REGION_CODING_0],
                    fivePrimePhaseCounts[GENE_PHASING_REGION_CODING_1], fivePrimePhaseCounts[GENE_PHASING_REGION_CODING_2],
                    threePrimePhaseCounts[GENE_PHASING_REGION_5P_UTR], threePrimePhaseCounts[GENE_PHASING_REGION_CODING_0],
                    threePrimePhaseCounts[GENE_PHASING_REGION_CODING_1], threePrimePhaseCounts[GENE_PHASING_REGION_CODING_2]);
        }
    }



    public static int phaseToRegion(int phase)
    {
        switch(phase)
        {
            case -1: return GENE_PHASING_REGION_5P_UTR;
            case 0: return GENE_PHASING_REGION_CODING_0;
            case 1: return GENE_PHASING_REGION_CODING_1;
            case 2: return GENE_PHASING_REGION_CODING_2;
        }

        return GENE_PHASING_REGION_5P_UTR;
    }

    public static int regionToPhase(int region)
    {
        if(region == GENE_PHASING_REGION_5P_UTR) return -1;
        if(region == GENE_PHASING_REGION_CODING_0) return 0;
        if(region == GENE_PHASING_REGION_CODING_1) return 1;
        if(region == GENE_PHASING_REGION_CODING_2) return 2;
        if(region == GENE_PHASING_REGION_PROMOTOR) return -1;

        return -1;
    }

    public void writeGeneLikelihoodData(final String outputDir)
    {
        try
        {
            String outputFilename = outputDir + "GFL_GENE_PHASE_COUNTS.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("GeneId,GeneName,Chromosome,Arm,GeneStart,GeneEnd,Strand");
            writer.write(",FivePrimeUTR,FivePrimePhase0,FivePrimePhase1,FivePrimePhase2");
            writer.write(",ThreePrimeUTR,ThreePrimePhase0,ThreePrimePhase1,ThreePrimePhase2");
            writer.newLine();

            for(Map.Entry<String, List<GeneRangeData>> entry : mChrForwardGeneDataMap.entrySet())
            {
                for(final GeneRangeData geneData :entry.getValue())
                {
                    if(!geneData.hasCodingTranscripts(true))
                        continue;

                    writer.write(String.format("%s,%s,%s,%s,%d,%d,%d",
                            geneData.GeneData.GeneId, geneData.GeneData.GeneName, geneData.GeneData.Chromosome, geneData.Arm,
                            geneData.GeneData.GeneStart, geneData.GeneData.GeneEnd, geneData.GeneData.Strand));

                    for(int i = 0; i <= 1; ++i)
                    {
                        final int[] phaseCounts = geneData.getPhaseCounts(i == 0);

                        writer.write(String.format(",%d,%d,%d,%d",
                                phaseCounts[GENE_PHASING_REGION_5P_UTR], phaseCounts[GENE_PHASING_REGION_CODING_0],
                                phaseCounts[GENE_PHASING_REGION_CODING_1], phaseCounts[GENE_PHASING_REGION_CODING_2]));

                        writer.newLine();
                    }
                }
            }

            closeBufferedWriter(writer);
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing disruptions: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addCmdLineArgs(options);
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Ensembl gene transcript data cache directory");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        Configurator.setRootLevel(Level.DEBUG);

        LOGGER.info("Generating gene likelihood data");

        SvGeneTranscriptCollection ensemblDataCache = new SvGeneTranscriptCollection();
        ensemblDataCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));

        List<String> restrictedGeneIds = Lists.newArrayList();
        if(cmd.hasOption(LIMITED_GENE_IDS))
        {
            restrictedGeneIds = Arrays.stream(cmd.getOptionValue(LIMITED_GENE_IDS).split(";")).collect(Collectors.toList());
        }

        boolean limitedLoading = !restrictedGeneIds.isEmpty();

        if(!ensemblDataCache.loadEnsemblData(limitedLoading))
        {
            LOGGER.error("Ensembl data cache load failed, exiting");
            return;
        }

        if(limitedLoading)
        {
            ensemblDataCache.loadEnsemblTranscriptData(restrictedGeneIds);
        }

        String outputDir = formOutputPath(cmd.getOptionValue(DATA_OUTPUT_DIR));
        FusionLikelihood fusionLikelihood = new FusionLikelihood();
        fusionLikelihood.initialise(cmd, ensemblDataCache);
        fusionLikelihood.generateGenePhasingCounts();
        fusionLikelihood.writeGeneLikelihoodData(outputDir);

        LOGGER.info("Gene likelihood data generation complete");
    }

    /*
        public void investigateExonOverlaps(final EnsemblGeneData geneData, final List<TranscriptExonData> transExonDataList)
    {
        int exonMatchDiffPhases = 0;
        int exonMatches = 0;
        int overlapsSamePhases = 0;
        int overlapsDiffPhases = 0;
        List<Integer> uniqueTrans = Lists.newArrayList();

        for(int i = 0; i < transExonDataList.size(); ++i)
        {
            final TranscriptExonData exonData1 = transExonDataList.get(i);

            if(!uniqueTrans.contains(exonData1.TransId))
                uniqueTrans.add(exonData1.TransId);

            for(int j = i+1; j < transExonDataList.size(); ++j)
            {
                final TranscriptExonData exonData2 = transExonDataList.get(j);

                if(exonData2.TransId == exonData1.TransId)
                    continue;

                if(exonData1.ExonEnd < exonData2.ExonStart || exonData1.ExonStart > exonData2.ExonEnd)
                    continue;

                if(exonData1.ExonPhase == -1 && exonData2.ExonPhase == -1 && exonData1.ExonPhaseEnd == -1 && exonData2.ExonPhaseEnd == -1)
                    continue;

                    // there's an overlap - check phasing and exact base matches
                if(exonData1.ExonStart == exonData2.ExonStart && exonData1.ExonEnd == exonData2.ExonEnd)
                {
                    if(exonData1.ExonPhase == exonData2.ExonPhase && exonData1.ExonPhaseEnd == exonData2.ExonPhaseEnd)
                    {
                        ++exonMatches;
                    }
                    else
                    {
                        ++exonMatchDiffPhases;
                    }
                }
                else
                {
                    if(exonData1.ExonPhase == exonData2.ExonPhase && exonData1.ExonPhaseEnd == exonData2.ExonPhaseEnd)
                    {
                        ++overlapsSamePhases;
                    }
                    else
                    {
                        ++overlapsDiffPhases;
                    }
                }
            }
        }

        LOGGER.debug("gene({}: {}) transcripts({}) exons({}) match(samePhase={} diffPhase={}) overlap(samePhase={} diffPhase={})",
                geneData.GeneId, geneData.GeneName, uniqueTrans.size(), transExonDataList.size(),
                exonMatches, exonMatchDiffPhases, overlapsSamePhases, overlapsDiffPhases);
    }
     */

}
