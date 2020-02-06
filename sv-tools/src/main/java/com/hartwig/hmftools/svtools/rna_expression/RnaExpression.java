package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TRANS_COUNT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.UNIQUE_TRANS_COUNT;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.SAMPLE;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.createCmdLineOptions;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
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

public class RnaExpression
{
    private final RnaExpConfig mConfig;
    private final String mSampledId;
    private final RnaBamReader mRnaBamReader;
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final GcBiasAdjuster mGcBiasAdjuster;
    private final Map<Integer,Integer> mFragmentLengths;

    private BufferedWriter mWriter;
    private BufferedWriter mExonDataWriter;

    private static final Logger LOGGER = LogManager.getLogger(RnaExpression.class);

    public RnaExpression(final CommandLine cmd)
    {
        mConfig = new RnaExpConfig(cmd);

        mRnaBamReader = new RnaBamReader(mConfig);
        mGcBiasAdjuster = new GcBiasAdjuster(mConfig);

        mSampledId = cmd.getOptionValue(SAMPLE);

        mGeneTransCache = new SvGeneTranscriptCollection();
        mGeneTransCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));

        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            mGeneTransCache.setRestrictedGeneIdList(mConfig.RestrictedGeneIds);
        }

        mGeneTransCache.setRequiredData(true, false, false, !mConfig.AllTranscripts);
        mGeneTransCache.loadEnsemblData(false);

        mFragmentLengths = Maps.newHashMap();

        mWriter = null;
        mExonDataWriter = null;
    }

    public void runAnalysis()
    {
        if(!mRnaBamReader.validReader())
        {
            LOGGER.warn("BAM reader init failed");
            return;
        }

        if(mGcBiasAdjuster.enabled())
        {
            mGcBiasAdjuster.loadData();
            mGcBiasAdjuster.generateDepthCounts(mRnaBamReader, mGeneTransCache.getChrGeneDataMap());
        }

        // measure read counts of exonic regions for all specific genes
        int geneCount = 0;
        for(Map.Entry<String,List<EnsemblGeneData>> entry : mGeneTransCache.getChrGeneDataMap().entrySet())
        {
            final List<EnsemblGeneData> genesDataList = entry.getValue();

            for(EnsemblGeneData geneData : genesDataList)
            {
                processGene(geneData);
                ++geneCount;

                if(geneCount > 1 && (geneCount % 100) == 0)
                    LOGGER.info("processed {} genes", geneCount);
            }
        }

        if(mConfig.WriteFragmentLengths && !mFragmentLengths.isEmpty())
        {
            LOGGER.info("fragment sizes:");

            for (Map.Entry<Integer, Integer> entry : mFragmentLengths.entrySet())
            {
                LOGGER.info("FRAG_SIZE: {},{}", entry.getKey(), entry.getValue());
            }
        }

        closeBufferedWriter(mWriter);
        closeBufferedWriter(mExonDataWriter);
        mRnaBamReader.close();
    }

    private void processGene(final EnsemblGeneData geneData)
    {
        GeneReadData geneReadData = new GeneReadData(geneData);

        List<TranscriptData> transDataList = Lists.newArrayList(mGeneTransCache.getTranscripts(geneData.GeneId));

        if(transDataList.isEmpty())
        {
            LOGGER.warn("no transcripts found for gene({}:{})", geneData.GeneId, geneData.GeneName);
            return;
        }

        if(!mConfig.SpecificTransIds.isEmpty())
            transDataList = transDataList.stream().filter(x -> mConfig.SpecificTransIds.contains(x.TransName)).collect(Collectors.toList());

        // form a genomic region for each unique exon amongst the transcripts
        long minTransPos = -1;
        long maxTransPos = 0;

        for(final TranscriptData transData : transDataList)
        {
            RegionReadData prevRegionReadData = null;

            for(int i = 0; i < transData.exons().size(); ++ i)
            {
                ExonData exon = transData.exons().get(i);

                RegionReadData exonReadData = geneReadData.findExonRegion(exon.ExonStart, exon.ExonEnd);

                if (exonReadData == null)
                {
                    GenomeRegion region = GenomeRegions.create(geneData.Chromosome, exon.ExonStart, exon.ExonEnd);
                    exonReadData = new RegionReadData(region);
                    geneReadData.addExonRegion(exonReadData);
                }

                exonReadData.addExonRef(transData.TransName, exon.ExonRank);

                if(prevRegionReadData != null)
                {
                    prevRegionReadData.addPostRegion(exonReadData);
                    exonReadData.addPreRegion(prevRegionReadData);
                }

                // create intronic regions
                if(prevRegionReadData != null)
                {
                    long intronStart = prevRegionReadData.end() + 1;
                    long intronEnd = exon.ExonStart - 1;

                    RegionReadData intronReadData = geneReadData.createOrFindIntronRegion(intronStart, intronEnd);
                    intronReadData.addExonRef(transData.TransName, exon.ExonRank);
                }

                prevRegionReadData = exonReadData;
            }

            maxTransPos = max(transData.TransEnd, maxTransPos);

            if(minTransPos < 0 || transData.TransStart < minTransPos)
                minTransPos = transData.TransStart;
        }

        // cache reference bases for comparison with read bases
        if(mConfig.RefFastaSeqFile != null)
        {
            for (RegionReadData region : geneReadData.getExonRegions())
            {
                final String regionRefBases = mConfig.RefFastaSeqFile.getSubsequenceAt(
                        region.chromosome(), region.start(), region.end()).getBaseString();

                region.setRefBases(regionRefBases);
            }
        }

        // use a buffer around the gene to pick up reads which span outside its transcripts
        GenomeRegion geneRegion = GenomeRegions.create(
                geneData.Chromosome, minTransPos - mConfig.MaxFragmentSize, maxTransPos + mConfig.MaxFragmentSize);

        mRnaBamReader.readBamCounts(geneReadData, geneRegion);
        mRnaBamReader.analyseReads();

        // report evidence for each gene transcript
        for(final TranscriptData transData : transDataList)
        {
            final TranscriptResults results = calculateTranscriptResults(geneReadData, transData);
            geneReadData.getTranscriptResults().add(results);
            writeResults(geneReadData, results);

            if(mConfig.WriteExonData)
            {
                writeExonData(geneReadData, transData);
            }

            for(Integer fragmentLength : geneReadData.getFragmentLengths())
            {
                fragmentLength = min(fragmentLength, 5000); // to prevent map blowing out in size

                Integer count = mFragmentLengths.get(fragmentLength);
                if(count == null)
                    mFragmentLengths.put(fragmentLength, 1);
                else
                    mFragmentLengths.put(fragmentLength, count + 1);
            }
        }

        // mGeneReadDatalist.add(geneReadData);
    }

    private TranscriptResults calculateTranscriptResults(final GeneReadData geneReadData, final TranscriptData transData)
    {
        int exonsFound = 0;

        int fragments = 0;
        int uniqueFragments = 0;

        int sjFragments = 0;
        int sjUniqueFragments = 0;
        int spliceJunctionsSupported = 0;

        long exonicBases = 0;
        long exonicBaseCoverage = 0;

        /* Criteria for transcript selection
        - all exon junctions covered
        - unique exon junctions
        - split reads skipping exons
        - unique exon reads (but could cover introns as well
         */

        final List<ExonData> exons = transData.exons();

        for(int i = 0; i < exons.size(); ++i)
        {
            ExonData exon = exons.get(i);

            final RegionReadData exonReadData = geneReadData.findExonRegion(exon.ExonStart, exon.ExonEnd);
            if(exonReadData == null)
                continue;

            int exonCoverage = exonReadData.baseCoverage(1);
            exonicBaseCoverage += exonCoverage;

            if(exonCoverage > 0)
                ++exonsFound;

            exonicBases += exon.ExonEnd - exon.ExonStart + 1;

            boolean linked = true;

            if(i > 0)
            {
                int[] sjReads = exonReadData.getTranscriptJunctionMatchCount(transData.TransName, SE_START);

                if(sjReads[TRANS_COUNT] == 0)
                {
                    linked = false;
                }
                else
                {
                    sjFragments += sjReads[TRANS_COUNT];
                    sjUniqueFragments += sjReads[UNIQUE_TRANS_COUNT];
                }
            }

            if(linked)
                ++spliceJunctionsSupported;
        }

        int[] supportingFragments = geneReadData.getTranscriptReadCount(transData.TransName);

        TranscriptResults results = ImmutableTranscriptResults.builder()
                .trans(transData)
                .exonsFound(exonsFound)
                .supportingFragments(supportingFragments[TRANS_COUNT])
                .uniqueFragments(supportingFragments[UNIQUE_TRANS_COUNT])
                .spliceJunctionsSupported(spliceJunctionsSupported)
                .spliceJunctionFragments(sjFragments)
                .spliceJunctionUniqueFragments(sjUniqueFragments)
                .exonicBases(exonicBases)
                .exonicBaseCoverage(exonicBaseCoverage)
                .build();

        return results;
    }

    private void writeResults(final GeneReadData geneReadData, final TranscriptResults transResults)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mWriter == null)
            {
                final String outputFileName = mConfig.OutputDir + "RNA_GENE_EXPRESSION.csv";

                mWriter = createBufferedWriter(outputFileName, false);
                mWriter.write("SampleId,GeneId,GeneName,TransId,Canonical,ExonCount");
                mWriter.write(",ExonsMatched,ExonicBases,ExonicCoverage,Fragments,UniqueFragments");
                mWriter.write(",SpliceJuncSupported,SpliceJuncFragments,UniqueSpliceJuncFragments");
                mWriter.newLine();
            }

            final TranscriptData transData = transResults.trans();

            mWriter.write(String.format("%s,%s,%s,%s,%s,%d",
                    mSampledId, geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName,
                    transData.TransName, transData.IsCanonical, transData.exons().size()));

            mWriter.write(String.format(",%d,%d,%d",
                    transResults.exonsFound(), transResults.exonicBases(), transResults.exonicBaseCoverage()));

            mWriter.write(String.format(",%d,%d,%d,%d,%d",
                    transResults.supportingFragments(), transResults.uniqueFragments(), transResults.spliceJunctionsSupported(),
                    transResults.spliceJunctionFragments(), transResults.spliceJunctionUniqueFragments()));

            mWriter.newLine();

        }
        catch(IOException e)
        {
            LOGGER.error("failed to write gene expression file: {}", e.toString());
        }
    }

    private void writeExonData(final GeneReadData geneReadData, final TranscriptData transData)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mExonDataWriter == null)
            {
                final String outputFileName = mConfig.OutputDir + "RNA_EXON_EXPRESSION.csv";

                mExonDataWriter = createBufferedWriter(outputFileName, false);
                mExonDataWriter.write("SampleId,GeneId,GeneName,TransId,ExonRank,ExonStart,ExonEnd");
                mExonDataWriter.write(",TotalCoverage,AvgDepth,Fragments,UniqueFragments");
                mExonDataWriter.write(",SpliceJuncStart,SpliceJuncEnd,UniqueSpliceJuncStart,UniqueSpliceJuncEnd");
                mExonDataWriter.newLine();
            }

            final List<ExonData> exons = transData.exons();

            for(int i = 0; i < exons.size(); ++i)
            {
                ExonData exon = exons.get(i);

                final RegionReadData exonReadData = geneReadData.findExonRegion(exon.ExonStart, exon.ExonEnd);
                if (exonReadData == null)
                    continue;

                mExonDataWriter.write(String.format("%s,%s,%s,%s,%d,%d,%d",
                        mSampledId, geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName,
                        transData.TransName, exon.ExonRank, exon.ExonStart, exon.ExonEnd));

                int[] matchCounts = exonReadData.getTranscriptReadCount(transData.TransName);
                int[] startSjCounts = exonReadData.getTranscriptJunctionMatchCount(transData.TransName, SE_START);
                int[] endSjCounts = exonReadData.getTranscriptJunctionMatchCount(transData.TransName, SE_END);

                mExonDataWriter.write(String.format(",%d,%.0f,%d,%d,%d,%d,%d,%d",
                        exonReadData.baseCoverage(1), exonReadData.averageDepth(),
                        matchCounts[TRANS_COUNT], matchCounts[UNIQUE_TRANS_COUNT],
                        startSjCounts[TRANS_COUNT], endSjCounts[TRANS_COUNT],
                        startSjCounts[UNIQUE_TRANS_COUNT], endSjCounts[UNIQUE_TRANS_COUNT]));

                /* match types are not recorded per transcript
                mExonDataWriter.write(String.format(",%d,%.0f,%d,%d,%d,%d,%d",
                        exonReadData.matchedReadCount(MATCH_TYPE_WITHIN_EXON),
                        exonReadData.matchedReadCount(MATCH_TYPE_EXON_BOUNDARY) + exonReadData.matchedReadCount(MATCH_TYPE_EXON_MATCH),
                        exonReadData.matchedReadCount(MATCH_TYPE_UNSPLICED));
                 */

                mExonDataWriter.newLine();
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write exon expression file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createCmdLineOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        if(!RnaExpConfig.checkValid(cmd))
        {
            LOGGER.error("missing config options, exiting");
            return;
        }

        RnaExpression rnaExpression = new RnaExpression(cmd);
        rnaExpression.runAnalysis();

        LOGGER.info("RNA expression analysis complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
