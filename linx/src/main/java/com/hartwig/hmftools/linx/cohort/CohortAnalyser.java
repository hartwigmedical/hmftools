package com.hartwig.hmftools.linx.cohort;

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
import static com.hartwig.hmftools.linx.cohort.SvType.classifySv;
import static com.hartwig.hmftools.linx.cohort.SvType.ignoreResolvedType;
import static com.hartwig.hmftools.linx.types.SvVarData.GENE_DATA_ITEM_DELIM;
import static com.hartwig.hmftools.linx.types.SvVarData.SV_DISRUPTIVE_STR;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.utils.StartEndPair;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.linx.types.ResolvedType;

import org.jetbrains.annotations.NotNull;

public class CohortAnalyser
{
    public final CohortConfig mConfig;

    // private final Map<String, List<SvData>> mSampleSvMap;

    private final PonCache mPonCache;
    private final boolean mBuildPon;

    // output writers
    private final BufferedWriter mSvDataWriter;
    private final BufferedWriter mSampleDataWriter;
    private final BufferedWriter mBreakendDataWriter;

    public CohortAnalyser(final ConfigBuilder configBuilder)
    {
        mConfig = new CohortConfig(configBuilder);

        // mSampleSvMap = Maps.newHashMap();
        mSvDataWriter = initialiseSvDataWriter();
        mSampleDataWriter = initialiseSampleDataWriter();
        mBreakendDataWriter = initialiseBreakendWriter();

        mPonCache = new PonCache();

        if(mConfig.PonFile != null)
        {
            mPonCache.loadPonFile(mConfig.PonFile);
            mBuildPon = false;
        }
        else
        {
            mBuildPon = mConfig.WriteTypes.contains(WriteType.PON);
        }
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

        if(mBuildPon)
            mPonCache.writePonCache(mConfig);

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
            // int Index = fieldsIndexMap.get(FLD_);
            // int Index = fieldsIndexMap.get(FLD_);

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
                    // mSampleSvMap.put(sampleId, svDataList);
                }

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

                svDataList.add(svData);

                svData.setSvType(classifySv(svData));
                svData.setGeneDisruptive();

                if(mBuildPon)
                {
                    if(svData.SvType == SGL)
                    {
                        mPonCache.registerSgl(svData.ChrStart, svData.PosStart, svData.OrientStart);
                    }
                    else
                    {
                        mPonCache.registerSv(
                                svData.ChrStart, svData.PosStart, svData.OrientStart,
                                svData.ChrEnd, svData.PosEnd, svData.OrientEnd);
                    }
                }
                else if(mPonCache.hasEntries())
                {
                    if(isSgl)
                    {
                        svData.setPonMatch(mPonCache.matchesSglEntry(
                                svData.ChrStart, svData.PosStart, svData.OrientStart, mConfig.PonMargin));
                    }
                    else
                    {
                        svData.setPonMatch(mPonCache.matchesSvEntry(
                                svData.ChrStart, svData.PosStart, svData.OrientStart, svData.ChrEnd, svData.PosEnd, svData.OrientEnd,
                                mConfig.PonMargin));
                    }
                }

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

    private BufferedWriter initialiseSvDataWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.SAMPLE_SUMMARY))
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
            sj.add(FLD_GENE_START).add(FLD_GENE_END);

            // sj.add(FLD_HOMOLOGY_START).add(FLD_HOMOLOGY_END).add(FLD_INS_SEQ);
            // sj.add("QualScore").add(FLD_AF_START).add(FLD_AF_END);
            sj.add(FLD_REPEAT_CLASS).add(FLD_PON_MATCH);

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

    private void writeTrimmedSvData(final String sampleId, final SvData svData) throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(sampleId).add(String.valueOf(svData.Id)).add(svData.SvType.toString());
        sj.add(svData.ChrStart).add(String.valueOf(svData.PosStart)).add(svData.OrientStart.toString());
        sj.add(svData.ChrEnd).add(String.valueOf(svData.PosEnd)).add(svData.SvType == SGL ? "0" : svData.OrientEnd.toString());
        sj.add(String.valueOf(svData.ClusterId));
        sj.add(String.valueOf(svData.ClusterCount));
        sj.add(String.valueOf(svData.ChainId));
        sj.add(String.valueOf(svData.ChainCount));
        sj.add(svData.ChainIndex);
        sj.add(svData.GeneStart);
        sj.add(svData.GeneEnd);
        sj.add(svData.RepeatClass);
        sj.add(svData.ponMatch().toString());

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

            for(SvType svType : SvType.values())
            {
                sj.add(svType.toString());
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

        for(SvType svType : SvType.values())
        {
            sj.add(String.valueOf(statistics.TypeCounts[svType.ordinal()]));
        }

        sj.add(String.valueOf(statistics.MeiLINE));
        sj.add(String.valueOf(statistics.MeiSINE));
        sj.add(String.valueOf(statistics.MeiOther));
        sj.add(String.valueOf(statistics.Genic));
        sj.add(String.valueOf(statistics.GeneDisruptive));

        mSampleDataWriter.write(sj.toString());
        mSampleDataWriter.newLine();
    }

    private BufferedWriter initialiseBreakendWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.SAMPLE_SUMMARY))
            return null;

        try
        {
            String filename = mConfig.formOutputFilename("breakend");
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            // definitional fields
            sj.add(FLD_SAMPLE_ID).add(FLD_SV_ID).add(FLD_TYPE);
            sj.add(FLD_CHROMOSOME).add(FLD_POSITION).add(FLD_ORIENTATION).add(FLD_PON_MATCH);
            sj.add("IsStart").add("PairedGene");

            sj.add(FLD_CLUSTER_COUNT).add(FLD_CHAIN_COUNT);

            sj.add(FLD_GENE_NAME).add(FLD_REGION_TYPE).add(FLD_CODING_TYPE);
            sj.add("Disruptive").add("Exon");

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
        if(mBreakendDataWriter == null)
            return;

        if(svData.GeneStart.isEmpty() && svData.GeneEnd.isEmpty())
            return;

        boolean isPaired = !svData.GeneStart.isEmpty() && !svData.GeneEnd.isEmpty();

        for(StartEndPair se : StartEndPair.values())
        {
            String genesData = se.isStart() ? svData.GeneStart : svData.GeneEnd;

            if((se.isEnd() && svData.SvType == SGL) || genesData.isEmpty())
                continue;

            String[] geneDataItems = genesData.split(ITEM_DELIM, -1);

            for(String geneData : geneDataItems)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(sampleId).add(String.valueOf(svData.Id)).add(svData.SvType.toString());

                if(se.isStart())
                    sj.add(svData.ChrStart).add(String.valueOf(svData.PosStart)).add(svData.OrientStart.toString());
                else
                    sj.add(svData.ChrEnd).add(String.valueOf(svData.PosEnd)).add(svData.OrientEnd.toString());

                sj.add(svData.ponMatch().toString());
                sj.add(String.valueOf(se.isStart()));
                sj.add(String.valueOf(isPaired));

                sj.add(String.valueOf(svData.ClusterCount));
                sj.add(String.valueOf(svData.ChainCount));

                // format: geneName|regionType|codingType|disruption or not present|exon=123 or not present
                String[] geneParts = geneData.split("\\" + GENE_DATA_ITEM_DELIM, -1);

                if(geneParts.length < 3)
                    continue;

                int partIndex = 0;
                String geneName = geneParts[partIndex++];
                String regionType = geneParts[partIndex++];
                String codingType = geneParts[partIndex++];

                sj.add(geneName);
                sj.add(regionType);
                sj.add(codingType);

                boolean isDisruptive = false;

                if(partIndex < geneParts.length && geneParts[partIndex].equals(SV_DISRUPTIVE_STR))
                {
                    ++partIndex;
                    isDisruptive = true;
                }

                sj.add(String.valueOf(isDisruptive));

                int exon = -1;
                if(partIndex < geneParts.length && geneParts[partIndex].contains("exon"))
                {
                    String[] exonData = geneParts[partIndex].split("=", 2);
                    exon = Integer.parseInt(exonData[1]);
                }

                sj.add(String.valueOf(exon));

                mBreakendDataWriter.write(sj.toString());
                mBreakendDataWriter.newLine();
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
