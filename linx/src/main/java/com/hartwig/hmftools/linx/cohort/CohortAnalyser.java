package com.hartwig.hmftools.linx.cohort;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHR_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHR_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CODING_TYPE;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENT_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENT_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_TYPE;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TYPE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxApplication.APP_NAME;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.COHORT_FILE_ID_GERMLINE_SVS;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_AF_END;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_AF_START;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_CHAIN_COUNT;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_CHAIN_ID;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_CHAIN_INDEX;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_CLUSTER_COUNT;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_CLUSTER_ID;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_GENE_END;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_GENE_START;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_LINE_END;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_LINE_START;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_REPEAT_CLASS;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_RESOLVED_TYPE;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.FLD_SV_ID;
import static com.hartwig.hmftools.linx.cohort.CohortDataWriter.cohortDataFilename;
import static com.hartwig.hmftools.linx.cohort.PonMatchType.NONE;
import static com.hartwig.hmftools.linx.cohort.SvCategory.classifySv;
import static com.hartwig.hmftools.linx.cohort.SvCategory.ignoreResolvedType;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.utils.StartEndPair;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.linx.types.ResolvedType;

import org.jetbrains.annotations.NotNull;

public class CohortAnalyser
{
    public final CohortConfig mConfig;

    private final PonCache mPonCache;
    private final Map<String,PonCache> mPonCacheMap;
    private final PonBuilder mPonBuilder;

    // output writers
    private final BufferedWriter mSvDataWriter;
    private final BufferedWriter mSampleDataWriter;
    private final BufferedWriter mBreakendDataWriter;
    private final BufferedWriter mDisruptiveBreakendDataWriter;

    public CohortAnalyser(final ConfigBuilder configBuilder)
    {
        mConfig = new CohortConfig(configBuilder);

        mPonCacheMap = Maps.newHashMap();

        if(mConfig.PonFile != null)
        {
            if(mConfig.PonFile.contains(ITEM_DELIM))
            {
                String[] ponFiles = mConfig.PonFile.split(ITEM_DELIM, -1);

                for(String ponFileEntry : ponFiles)
                {
                    String[] parts = ponFileEntry.split("=", 2);

                    if(parts.length != 2)
                    {
                        LNX_LOGGER.error("invalid PON file entry: {}", ponFileEntry);
                        System.exit(1);
                    }

                    String ponFilename = parts[1];
                    PonCache ponCache = new PonCache();
                    ponCache.loadPonFile(ponFilename);
                    mPonCacheMap.put(parts[0], ponCache);
                }

                mPonCache = null;
            }
            else
            {
                mPonCache = new PonCache();
                mPonCache.loadPonFile(mConfig.PonFile);
            }
        }
        else
        {
            mPonCache = null;
        }

        if(mConfig.WriteTypes.contains(WriteType.PON))
        {
            mPonBuilder = new PonBuilder();
        }
        else
        {
            mPonBuilder = null;
        }

        mSvDataWriter = initialiseSvDataWriter();
        mSampleDataWriter = initialiseSampleDataWriter();
        mBreakendDataWriter = initialiseBreakendWriter(false);
        mDisruptiveBreakendDataWriter = initialiseBreakendWriter(true);
    }

    public void run()
    {
        LNX_LOGGER.info("start cohort analysis");

        long startTimeMs = System.currentTimeMillis();

        String svCohortFile = cohortDataFilename(mConfig.InputDir, COHORT_FILE_ID_GERMLINE_SVS);
        loadSvData(svCohortFile);

        closeBufferedWriter(mSvDataWriter);
        closeBufferedWriter(mSampleDataWriter);
        closeBufferedWriter(mBreakendDataWriter);
        closeBufferedWriter(mDisruptiveBreakendDataWriter);

        if(mPonBuilder != null)
            mPonBuilder.writePonCache(mConfig);

        LNX_LOGGER.info("completed cohort analysis, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void loadSvData(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String header = fileReader.readLine();
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int sampleIndex = fieldsIndexMap.get(FLD_SAMPLE_ID);
            int svIdIndex = fieldsIndexMap.get(FLD_SV_ID);
            int chrStartIndex = fieldsIndexMap.get(FLD_CHR_START);
            int chrEndIndex = fieldsIndexMap.get(FLD_CHR_END);
            int posStartIndex = fieldsIndexMap.get(FLD_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_POS_END);
            int orientStartIndex = fieldsIndexMap.get(FLD_ORIENT_START);
            int orientEndIndex = fieldsIndexMap.get(FLD_ORIENT_END);

            int typeIndex = fieldsIndexMap.get(FLD_TYPE);
            int resolvedIndex = fieldsIndexMap.get(FLD_RESOLVED_TYPE);
            int geneStartIndex = fieldsIndexMap.get(FLD_GENE_START);
            int geneEndIndex = fieldsIndexMap.get(FLD_GENE_END);
            int clusterIdIndex = fieldsIndexMap.get(FLD_CLUSTER_ID);
            int clusterCountIndex = fieldsIndexMap.get(FLD_CLUSTER_COUNT);
            int chainIdIndex = fieldsIndexMap.get(FLD_CHAIN_ID);
            int chainCountIndex = fieldsIndexMap.get(FLD_CHAIN_COUNT);
            int chainIndexIndex = fieldsIndexMap.get(FLD_CHAIN_INDEX);
            int lineStartIndex = fieldsIndexMap.get(FLD_LINE_START);
            int lineEndIndex = fieldsIndexMap.get(FLD_LINE_END);
            int afStartIndex = fieldsIndexMap.get(FLD_AF_START);
            int afEndIndex = fieldsIndexMap.get(FLD_AF_END);

            int repeatClassIndex = fieldsIndexMap.get(FLD_REPEAT_CLASS);
            // spare: int Index = fieldsIndexMap.get(FLD_);

            Set<String> restrictedSampleIds = !mConfig.SampleIds.isEmpty() ? mConfig.SampleIds.stream().collect(Collectors.toSet()) : null;

            String line;
            String currentSample = "";
            List<SvData> svDataList = null;

            int svCount = 0;
            int sampleCount = 0;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                String sampleId = values[sampleIndex];

                if(!sampleId.equals(currentSample))
                {
                    ++sampleCount;

                    if((sampleCount % 100) == 0)
                    {
                        LNX_LOGGER.debug("processed {} samples", sampleCount);
                    }

                    if(svDataList != null)
                    {
                        processSampleSvs(currentSample, svDataList);
                    }

                    currentSample = sampleId;
                    svDataList = Lists.newArrayList();
                }

                if(restrictedSampleIds != null && !restrictedSampleIds.contains(sampleId))
                    continue;

                ++svCount;

                if((svCount % 1_000_000) == 0)
                {
                    LNX_LOGGER.debug("processed {} SVs", svCount);
                }

                ResolvedType resolvedType = ResolvedType.valueOf(values[resolvedIndex]);

                if(ignoreResolvedType(resolvedType))
                    continue;

                boolean isLine = !values[lineStartIndex].equals("NONE") || !values[lineEndIndex].equals("NONE");
                double afStart = Double.parseDouble(values[afStartIndex]);
                double afEnd = Double.parseDouble(values[afEndIndex]);
                double averageAf = (afStart + afEnd) * 0.5;

                StructuralVariantType svType = StructuralVariantType.valueOf(values[typeIndex]);
                boolean isSgl = svType == SGL;

                int chainId = Integer.parseInt(values[chainIdIndex]);
                int chainCount = Integer.parseInt(values[chainCountIndex]);

                SvData svData = new SvData(
                        Integer.parseInt(values[svIdIndex]),
                        values[chrStartIndex], Integer.parseInt(values[posStartIndex]), Orientation.fromByteStr(values[orientStartIndex]),
                        !isSgl ? values[chrEndIndex] : "", !isSgl ? Integer.parseInt(values[posEndIndex]) : 0,
                        !isSgl ? Orientation.fromByteStr(values[orientEndIndex]) : null,
                        svType, resolvedType,
                        values[geneStartIndex], values[geneEndIndex],
                        Integer.parseInt(values[clusterIdIndex]), Integer.parseInt(values[clusterCountIndex]), isLine,
                        chainId, chainCount, chainCount > 0 ? values[chainIndexIndex] : "",
                        averageAf, values[repeatClassIndex]);

                svData.setSvType(classifySv(svData));
                svData.setGeneDisruptive();

                if(mPonBuilder != null)
                {
                    if(svData.SvType == SGL)
                    {
                        mPonBuilder.registerSgl(svData.ChrStart, svData.PosStart, svData.OrientStart);
                    }
                    else
                    {
                        mPonBuilder.registerSv(
                                svData.ChrStart, svData.PosStart, svData.OrientStart,
                                svData.ChrEnd, svData.PosEnd, svData.OrientEnd);
                    }
                }
                else
                {
                    markPonEntries(svData);

                    if(mConfig.RestrictNonPon && svData.hasPonMatch())
                        continue;
                }

                svDataList.add(svData);

                writeTrimmedSvData(sampleId, svData);
                writeBreakendData(sampleId, svData);
            }

            if(svDataList != null)
                processSampleSvs(currentSample, svDataList);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("failed to read cohort SV file({})", filename, exception.toString());
        }
    }

    private void markPonEntries(final SvData svData)
    {
        if(mPonCache != null && mPonCache.hasEntries())
        {
            PonMatch ponMatch = PonMatch.NONE;

            if(svData.SvType == SGL)
            {
                ponMatch = mPonCache.matchesSglEntry(
                        svData.ChrStart, svData.PosStart, svData.OrientStart, mConfig.PonMargin);
            }
            else
            {
                ponMatch = mPonCache.matchesSvEntry(
                        svData.ChrStart, svData.PosStart, svData.OrientStart, svData.ChrEnd, svData.PosEnd, svData.OrientEnd,
                        mConfig.PonMargin);
            }

            if(ponMatch != PonMatch.NONE)
                svData.addPonMatch(ponMatch);

            return;
        }

        if(!mPonCacheMap.isEmpty())
        {
            for(Map.Entry<String,PonCache> entry : mPonCacheMap.entrySet())
            {
                String ponName = entry.getKey();
                PonCache ponCache = entry.getValue();
                PonMatch ponMatch = PonMatch.NONE;

                if(svData.SvType == SGL)
                {
                    ponMatch = ponCache.matchesSglEntry(
                            svData.ChrStart, svData.PosStart, svData.OrientStart, mConfig.PonMargin);
                }
                else
                {
                    ponMatch = ponCache.matchesSvEntry(
                            svData.ChrStart, svData.PosStart, svData.OrientStart, svData.ChrEnd, svData.PosEnd, svData.OrientEnd,
                            mConfig.PonMargin);
                }

                if(ponMatch != PonMatch.NONE)
                    svData.addPonMatch(ponName, ponMatch);
            }
        }
    }

    private void processSampleSvs(final String sampleId, final List<SvData> svDataList) throws IOException
    {
        if(mSampleDataWriter == null)
            return;

        SvStatistics statistics = new SvStatistics();

        for(SvData sv : svDataList)
        {
            statistics.addSvData(sv);
        }

        writeSampleData(sampleId, statistics);
    }

    private static final String FLD_PON_MATCH = "PonMatch";
    private static final String FLD_CATEGORY = "Category";

    private BufferedWriter initialiseSvDataWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.SV_DATA))
            return null;

        try
        {
            String filename = mConfig.formOutputFilename("sv_data");
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_SAMPLE_ID).add(FLD_SV_ID).add(FLD_TYPE);
            sj.add(FLD_CHR_START).add(FLD_POS_START).add(FLD_ORIENT_START);
            sj.add(FLD_CHR_END).add(FLD_POS_END).add(FLD_ORIENT_END);

            sj.add(FLD_CLUSTER_ID).add(FLD_CLUSTER_COUNT);
            sj.add(FLD_CHAIN_ID).add(FLD_CHAIN_COUNT).add(FLD_CHAIN_INDEX);
            sj.add(FLD_RESOLVED_TYPE).add(FLD_CATEGORY);
            sj.add(FLD_GENE_START).add(FLD_GENE_END);

            // sj.add(FLD_HOMOLOGY_START).add(FLD_HOMOLOGY_END).add(FLD_INS_SEQ);
            // sj.add("QualScore").add(FLD_AF_START).add(FLD_AF_END);
            sj.add(FLD_REPEAT_CLASS);
            addPonColumns(sj);

            writer.write(sj.toString());

            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to write output file: {}", e.toString());
            return null;
        }
    }

    private void addPonColumns(final StringJoiner sj)
    {
        if(mPonCache != null)
        {
            sj.add(FLD_PON_MATCH);
        }
        else
        {
            for(String ponName : mPonCacheMap.keySet())
            {
                sj.add(format("%s_%s", FLD_PON_MATCH, ponName));
            }
        }
    }

    private void addPonValues(final SvData svData, final StringJoiner sj)
    {
        if(mPonCache != null)
        {
            PonMatch ponMatch = svData.getUnnamedPonMatch();
            addPonMatch(sj, ponMatch);
        }
        else
        {
            for(String ponName : mPonCacheMap.keySet())
            {
                PonMatch ponMatch = svData.ponMatches().get(ponName);
                addPonMatch(sj, ponMatch);
            }
        }
    }

    private void addPonMatch(final StringJoiner sj, final PonMatch ponMatch)
    {
        if(ponMatch == null || ponMatch == PonMatch.NONE)
        {
            if(mConfig.WritePonType)
                sj.add(NONE.toString());
            else
                sj.add("0");
        }
        else
        {
            if(mConfig.WritePonType)
                sj.add(ponMatch.toString());
            else
                sj.add(String.valueOf(ponMatch.Count));
        }
    }

    private void writeTrimmedSvData(final String sampleId, final SvData svData) throws IOException
    {
        if(mSvDataWriter == null)
            return;

        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(sampleId).add(String.valueOf(svData.Id)).add(svData.SvType.toString());
        sj.add(svData.ChrStart).add(String.valueOf(svData.PosStart)).add(svData.OrientStart.toString());
        sj.add(svData.ChrEnd).add(String.valueOf(svData.PosEnd)).add(svData.SvType == SGL ? "0" : svData.OrientEnd.toString());
        sj.add(String.valueOf(svData.ClusterId));
        sj.add(String.valueOf(svData.ClusterCount));
        sj.add(String.valueOf(svData.ChainId));
        sj.add(String.valueOf(svData.ChainCount));
        sj.add(svData.ChainIndex);
        sj.add(svData.Resolved.toString());
        sj.add(svData.category().toString());
        sj.add(svData.GeneStart);
        sj.add(svData.GeneEnd);
        sj.add(svData.RepeatClass);
        addPonValues(svData, sj);

        mSvDataWriter.write(sj.toString());
        mSvDataWriter.newLine();
    }

    private BufferedWriter initialiseSampleDataWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.SAMPLE_SUMMARY))
            return null;

        try
        {
            String filename = mConfig.formOutputFilename("sample_data");
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_SAMPLE_ID);
            sj.add("TotalCount");
            sj.add("NonPON");

            for(SvCategory svCategory : SvCategory.values())
            {
                sj.add(svCategory.toString());
            }

            sj.add("MeiLine");
            sj.add("MeiSine");
            sj.add("MeiOther");
            sj.add("Genic");
            sj.add("GeneDisruptive");

            writer.write(sj.toString());

            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to write output file: {}", e.toString());
            return null;
        }
    }

    private void writeSampleData(final String sampleId, final SvStatistics statistics) throws IOException
    {
        if(mSampleDataWriter == null)
            return;

        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(sampleId);
        sj.add(String.valueOf(statistics.Total));
        sj.add(String.valueOf(statistics.NonPON));

        for(SvCategory svCategory : SvCategory.values())
        {
            sj.add(String.valueOf(statistics.TypeCounts[svCategory.ordinal()]));
        }

        sj.add(String.valueOf(statistics.MeiLINE));
        sj.add(String.valueOf(statistics.MeiSINE));
        sj.add(String.valueOf(statistics.MeiOther));
        sj.add(String.valueOf(statistics.Genic));
        sj.add(String.valueOf(statistics.GeneDisruptive));

        mSampleDataWriter.write(sj.toString());
        mSampleDataWriter.newLine();
    }

    private BufferedWriter initialiseBreakendWriter(boolean disruptiveOnly)
    {
        if(disruptiveOnly && !mConfig.WriteTypes.contains(WriteType.DISRUPTIVE_BREAKEND))
            return null;

        if(!disruptiveOnly && !mConfig.WriteTypes.contains(WriteType.BREAKEND))
            return null;

        try
        {
            String fileId = disruptiveOnly ? "disruptive_breakend" : "breakend";
            String filename = mConfig.formOutputFilename(fileId);
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            // definitional fields
            sj.add(FLD_SAMPLE_ID).add(FLD_SV_ID).add(FLD_TYPE).add(FLD_CATEGORY);
            sj.add(FLD_CHROMOSOME).add(FLD_POSITION).add(FLD_ORIENTATION);

            sj.add("IsStart").add("PairedGene");

            sj.add(FLD_CLUSTER_COUNT).add(FLD_CHAIN_COUNT);

            sj.add(FLD_GENE_NAME).add(FLD_REGION_TYPE).add(FLD_CODING_TYPE);
            sj.add("Disruptive").add("Exon").add("Paired");

            addPonColumns(sj);

            writer.write(sj.toString());

            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to write output file: {}", e.toString());
            return null;
        }
    }

    private void writeBreakendData(final String sampleId, final SvData svData) throws IOException
    {
        if(mBreakendDataWriter == null && mDisruptiveBreakendDataWriter == null)
            return;

        if(svData.genesStart().isEmpty() && svData.genesEnd().isEmpty())
            return;

        for(StartEndPair se : StartEndPair.values())
        {
            if(se.isEnd() && svData.SvType == SGL)
                continue;

            List<GeneData> genesDataList = se.isStart() ? svData.genesStart() : svData.genesEnd();

            if(genesDataList.isEmpty())
                continue;

            for(GeneData geneData : genesDataList)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(sampleId).add(String.valueOf(svData.Id)).add(svData.SvType.toString()).add(svData.category().toString());

                if(se.isStart())
                    sj.add(svData.ChrStart).add(String.valueOf(svData.PosStart)).add(svData.OrientStart.toString());
                else
                    sj.add(svData.ChrEnd).add(String.valueOf(svData.PosEnd)).add(svData.OrientEnd.toString());

                sj.add(String.valueOf(se.isStart()));
                sj.add(String.valueOf(geneData.isPaired()));

                sj.add(String.valueOf(svData.ClusterCount));
                sj.add(String.valueOf(svData.ChainCount));

                sj.add(geneData.GeneName);
                sj.add(geneData.RegionType.toString());
                sj.add(geneData.CodingType.toString());

                boolean disruptive = geneData.Disruptive && !geneData.nonDisruptive();
                sj.add(String.valueOf(disruptive));
                sj.add(String.valueOf(geneData.Exon));
                sj.add(String.valueOf(geneData.isPaired()));

                addPonValues(svData, sj);

                if(mBreakendDataWriter != null)
                {
                    mBreakendDataWriter.write(sj.toString());
                    mBreakendDataWriter.newLine();
                }

                if(mDisruptiveBreakendDataWriter != null && disruptive)
                {
                    mDisruptiveBreakendDataWriter.write(sj.toString());
                    mDisruptiveBreakendDataWriter.newLine();
                }
            }
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CohortConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CohortAnalyser cohortAnalyser = new CohortAnalyser(configBuilder);
        cohortAnalyser.run();
    }
}
