package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.NONE;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.getHighestMatchType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class FusionWriter
{
    private final IsofoxConfig mConfig;
    private BufferedWriter mFusionWriter;
    private BufferedWriter mReadWriter;
    private BufferedWriter mFragmentWriter;
    private final boolean mWriteReads;

    public FusionWriter(final IsofoxConfig config)
    {
        mConfig = config;
        mWriteReads = mConfig.Fusions.WriteChimericReads;

        mFusionWriter = null;
        mReadWriter = null;
        mFragmentWriter = null;

        initialiseFusionWriter();
        initialiseReadWriter();
        initialiseFragmentWriter();
    }

    public void close()
    {
        closeBufferedWriter(mFusionWriter);
        closeBufferedWriter(mReadWriter);
        closeBufferedWriter(mFragmentWriter);
    }

    public static final String FUSION_FILE_ID = "fusions.csv";

    private void initialiseFusionWriter()
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile(FUSION_FILE_ID);

            mFusionWriter = createBufferedWriter(outputFileName, false);
            mFusionWriter.write(FusionReadData.csvHeader());
            mFusionWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to create fusions file: {}", e.toString());
        }
    }

    public synchronized void writeFusionData(Map<String,List<FusionReadData>> fusionCandidates)
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            for(Map.Entry<String, List<FusionReadData>> entry : fusionCandidates.entrySet())
            {
                for (final FusionReadData fusion : entry.getValue())
                {
                    mFusionWriter.write(fusion.toCsv());
                    mFusionWriter.newLine();
                }
            }

            if(mWriteReads)
            {
                for (List<FusionReadData> fusions : fusionCandidates.values())
                {
                    for (FusionReadData fusion : fusions)
                    {
                        for(Map.Entry<FusionFragmentType,List<FusionFragment>> entry : fusion.getFragments().entrySet())
                        {
                            for (FusionFragment fragment : entry.getValue())
                            {
                                writeFragmentData(fragment, fusionId(fusion.id()), entry.getKey());
                                writeReadData(fragment.reads(), fusionId(fusion.id()));
                            }
                        }
                    }
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusions file: {}", e.toString());
        }
    }

    private void initialiseReadWriter()
    {
        if(!mWriteReads)
            return;

        try
        {
            if(mReadWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("chimeric_reads.csv");

                mReadWriter = createBufferedWriter(outputFileName, false);
                mReadWriter.write("ReadSetCount,ReadId,FusionGroup,Chromosome,PosStart,PosEnd,Orientation,Cigar,InsertSize");
                mReadWriter.write(",FirstInPair,Supplementary,ReadReversed,ProperPair,SuppAlign");
                mReadWriter.write(",Bases,Flags,MateChr,MatePosStart,GeneSetStart,GeneSetEnd,TransExons,BestMatch,TransExonData");
                mReadWriter.newLine();
            }
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
            return;
        }
    }

    public void writeReadData(final List<ReadRecord> reads, final String groupStatus)
    {
        if(!mWriteReads)
            return;

        try
        {
            for(final ReadRecord read : reads)
            {
                mReadWriter.write(String.format("%s,%s,%s,%s,%d,%d,%d,%s,%d",
                        reads.size(), read.Id, groupStatus, read.Chromosome,
                        read.PosStart, read.PosEnd, read.orientation(), read.Cigar.toString(), read.fragmentInsertSize()));

                mReadWriter.write(String.format(",%s,%s,%s,%s,%s,%s,%d,%s,%d",
                        read.isFirstOfPair(), read.isSupplementaryAlignment(), read.isReadReversed(), read.isProperPair(),
                        read.hasSuppAlignment() ? read.getSuppAlignment().replaceAll(",", ";") : "NONE",
                        read.ReadBases, read.flags(), read.mateChromosome(), read.mateStartPosition()));

                // log the transcript exons affected, and the highest matching transcript
                String transExonData = "";

                int transExonRefCount = 0;
                RegionMatchType highestTransMatchType = getHighestMatchType(read.getTransExonRefs().keySet());

                if(highestTransMatchType != NONE)
                {
                    for (final TransExonRef transExonRef : read.getTransExonRefs().get(highestTransMatchType))
                    {
                        ++transExonRefCount;
                        transExonData = appendStr(transExonData, String.format("%s:%d:%s:%d",
                                transExonRef.GeneId, transExonRef.TransId, transExonRef.TransName, transExonRef.ExonRank), ';');
                    }
                }

                mReadWriter.write(String.format(",%d,%d,%d,%s,%s",
                        read.getGeneCollectons()[SE_START], read.getGeneCollectons()[SE_END],
                        transExonRefCount, highestTransMatchType, transExonRefCount == 0 ? "NONE" : transExonData));
                mReadWriter.newLine();
            }

        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
            return;
        }
    }

    public void writeUnfusedFragments(final Map<String,List<FusionFragment>> unfusedFragments)
    {
        if(!mWriteReads)
            return;

        for(List<FusionFragment> fragments : unfusedFragments.values())
        {
            for(FusionFragment fragment : fragments)
            {
                writeFragmentData(fragment, "UNFUSED", fragment.type());
                writeReadData(fragment.reads(), "UNFUSED");
            }
        }
    }

    private void initialiseFragmentWriter()
    {
        if(!mWriteReads)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile("chimeric_frags.csv");

            mFragmentWriter = createBufferedWriter(outputFileName, false);
            mFragmentWriter.write("ReadId,ReadCount,FusionGroup,Type,SameGene,ScCount");

            for(int se = SE_START; se <= SE_END; ++se)
            {
                final String prefix = se == SE_START ? "Start" : "End";
                mFragmentWriter.write(",Chr" + prefix);
                mFragmentWriter.write(",Orient" + prefix);
                mFragmentWriter.write(",JuncPos" + prefix);
                mFragmentWriter.write(",JuncOrient" + prefix);
                mFragmentWriter.write(",Region" + prefix);
                mFragmentWriter.write(",JuncType" + prefix);
            }

            mFragmentWriter.newLine();

        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric fragment data: {}", e.toString());
            return;
        }
    }

    public synchronized void writeFragmentData(final FusionFragment fragment, final String fusionId, FusionFragmentType type)
    {
        if(!mWriteReads)
            return;

        try
        {
            mFragmentWriter.write(String.format("%s,%d,%s,%s,%s,%d",
                    fragment.readId(), fragment.reads().size(), fusionId, type,
                    fragment.isSingleGene(), fragment.reads().stream().filter(x -> x.containsSoftClipping()).count()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                mFragmentWriter.write(String.format(",%s,%d,%d,%d,%s,%s",
                        fragment.chromosomes()[se], fragment.orientations()[se],
                        fragment.junctionPositions()[se], fragment.junctionOrientations()[se],
                        fragment.regionMatchTypes()[se], fragment.junctionTypes()[se]));
            }

            mFragmentWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric fragment data: {}", e.toString());
            return;
        }

    }

    public static Map<String,List<ReadRecord>> loadChimericReads(final String inputFile)
    {
        Map<String,List<ReadRecord>> readsMap = Maps.newHashMap();

        if(!Files.exists(Paths.get(inputFile)))
        {
            ISF_LOGGER.error("invalid chimeric reads file: {}", inputFile);
            return readsMap;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(inputFile));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                ISF_LOGGER.error("empty chimeric reads file: {}", inputFile);
                return readsMap;
            }

            final Map<String,Integer> fieldsMap  = createFieldsIndexMap(line, DELIMITER);
            String currentReadId = "";
            List<ReadRecord> currentReads = null;

            int readId = fieldsMap.get("ReadId");
            int chr = fieldsMap.get("Chromosome");
            int posStart = fieldsMap.get("PosStart");
            int posEnd = fieldsMap.get("PosEnd");
            int cigar = fieldsMap.get("Cigar");
            int insertSize = fieldsMap.get("InsertSize");
            int flags = fieldsMap.get("Flags");
            int suppAlgn = fieldsMap.get("SuppAlign");
            int bases = fieldsMap.get("Bases");
            int mateChr = fieldsMap.get("MateChr");
            int matePosStart = fieldsMap.get("MatePosStart");
            int geneSet = fieldsMap.get("GeneSet");
            int bestMatch = fieldsMap.get("BestMatch");
            int transExonData = fieldsMap.get("TransExonData");

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(DELIMITER, -1);

                ReadRecord read = new ReadRecord(
                        items[readId],
                        items[chr],
                        Integer.parseInt(items[posStart]),
                        Integer.parseInt(items[posEnd]),
                        items[bases],
                        cigarFromStr((items[cigar])),
                        Integer.parseInt(items[insertSize]),
                        Integer.parseInt(items[flags]),
                        items[mateChr],
                        Integer.parseInt(items[matePosStart]));

                String saData = items[suppAlgn];

                if(!saData.equals("NONE"))
                    read.setSuppAlignment(saData);

                read.setGeneCollection(SE_START,Integer.parseInt(items[geneSet]), true);
                read.setGeneCollection(SE_END,Integer.parseInt(items[geneSet]), true);

                RegionMatchType matchType = RegionMatchType.valueOf(items[bestMatch]);

                if(matchType != NONE)
                {
                    List<TransExonRef> transExonRefs = parseTransExonRefs(items[transExonData]);
                    read.getTransExonRefs().put(matchType, transExonRefs);
                }

                if(!read.Id.equals(currentReadId))
                {
                    currentReadId = read.Id;
                    currentReads = Lists.newArrayList();
                    readsMap.put(read.Id, currentReads);
                }

                currentReads.add(read);
            }

            ISF_LOGGER.info("loaded {} chimeric fragment reads from file({})", readsMap.size(), inputFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load chimeric reads file({}): {}", inputFile, e.toString());
        }

        return readsMap;
    }

    private static List<TransExonRef> parseTransExonRefs(final String data)
    {
        List<TransExonRef> transExonRefs = Lists.newArrayList();

        for(String ref : data.split(";", -1))
        {
            String[] items = ref.split(":");
            if(items.length != 4)
                continue;

            transExonRefs.add(new TransExonRef(items[0], Integer.parseInt(items[1]), items[2], Integer.parseInt(items[3])));
        }

        return transExonRefs;
    }

    private static Cigar cigarFromStr(final String cigarStr)
    {
        List<CigarElement> cigarElements = Lists.newArrayList();

        int index = 0;
        String basesStr = "";
        while(index < cigarStr.length())
        {
            char c = cigarStr.charAt(index);

            try
            {
                CigarOperator operator = CigarOperator.valueOf(String.valueOf(c));
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), operator));
                basesStr = "";
            }
            catch (Exception e)
            {
                basesStr += c;
            }

            /*
            if(c == CigarOperator.M.toString().charAt(0))
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), CigarOperator.M));
            else if(c == CigarOperator.I.toString().charAt(0))
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), CigarOperator.I));
            else if(c == CigarOperator.D.toString().charAt(0))
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), CigarOperator.D));
            else if(c == CigarOperator.N.toString().charAt(0))
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), CigarOperator.N));
            else if(c == CigarOperator.S.toString().charAt(0))
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), CigarOperator.S));
            else if(c == CigarOperator.H.toString().charAt(0))
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), CigarOperator.H));
            else if(c == CigarOperator.P.toString().charAt(0))
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), CigarOperator.P));
            else if(c == CigarOperator.X.toString().charAt(0))
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), CigarOperator.X));
            else
                basesStr += c;

            */

            ++index;
        }

        return new Cigar(cigarElements);
    }
}
