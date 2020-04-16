package com.hartwig.hmftools.isofox.novel;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionFinder
{
    private final IsofoxConfig mConfig;

    private final List<FusionReadData> mFusionReadData;
    private BufferedWriter mReadWriter;

    private final Map<String,List<ReadRecord>> mReadsMap;

    public FusionFinder(final IsofoxConfig config)
    {
        mConfig = config;

        mReadsMap = Maps.newHashMap();
        mFusionReadData = Lists.newArrayList();

        mReadWriter = null;
    }

    public void addChimericReads(final Map<String,List<ReadRecord>> chimericReadMap)
    {
        mergeChimericReadMaps(mReadsMap, chimericReadMap);
    }

    public static void addChimericRead(final Map<String,List<ReadRecord>> chimericReadMap, final ReadRecord read)
    {
        List<ReadRecord> chimericReads = chimericReadMap.get(read.Id);
        if (chimericReads == null)
        {
            chimericReads = Lists.newArrayList();
            chimericReadMap.put(read.Id, chimericReads);
        }

        chimericReads.add(read);
    }

    public static void mergeChimericReadMaps(final Map<String,List<ReadRecord>> destMap, final Map<String,List<ReadRecord>> sourceMap)
    {
        for(Map.Entry<String,List<ReadRecord>> entry :  sourceMap.entrySet())
        {
            List<ReadRecord> readsById = destMap.get(entry.getKey());

            if(readsById == null)
            {
                destMap.put(entry.getKey(), entry.getValue());
            }
            else
            {
                readsById.addAll(entry.getValue());
            }
        }
    }

    public void findFusions()
    {
        // make a list of chimeric reads by their ID
        int unpairedReads = 0;

        for(Map.Entry<String,List<ReadRecord>> entry : mReadsMap.entrySet())
        {
            final String readId = entry.getKey();

            final List<ReadRecord> reads = entry.getValue();

            if(reads.size() == 1)
            {
                ++unpairedReads;
                // continue;
            }

            // ISF_LOGGER.debug("read({}) count({})", readId, reads.size());

            writeReadData(reads);

            /*
            for(final ChimericRead read : reads)
            {
                ISF_LOGGER.debug("read({}) location({}:{} -> {}) cigar({}) dup({}) 2nd({}) supp({}) negStrand({})",
                        readId, read.Chromosome, read.PosStart, read.PosEnd, read.Cigar.toString(),
                        read.isDuplicate(), read.isSecondaryAlignment(), read.isSupplementaryAlignment(), read.isNegStrand());
            }
            */
        }

        ISF_LOGGER.info("chimeric fragments({}) unpaired({})", mReadsMap.size(), unpairedReads);


        closeBufferedWriter(mReadWriter);

        /*
        // first find read pairs if both reads are genic and of the correct orientations
        for(Map.Entry<String,List<ChimericRead>> entry : chimericReadsMap.entrySet())
        {
            final String chromosome = entry.getKey();

            final List<ChimericRead> chimericReads = entry.getValue();

            while(!chimericReads.isEmpty())
            {
                final ChimericRead read1 = chimericReads.get(0);
                final ChimericRead read2 =  findPairedRead(read1, chimericReadsMap);

                chimericReads.remove(0);

                if(read2 == null)
                    continue;




            }
        }

        */

    }

    public static final String FUSION_FILE_ID = "fusions.csv";

    private void writeFusionData()
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile(FUSION_FILE_ID);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write(AltSpliceJunction.csvHeader());
            writer.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusions file: {}", e.toString());
        }
    }

    private void writeReadData(final List<ReadRecord> reads)
    {
        try
        {
            if(mReadWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("chimeric_reads.csv");

                mReadWriter = createBufferedWriter(outputFileName, false);
                mReadWriter.write("ReadSetCount,ReadId,Chromosome,PosStart,PosEnd,Cigar,InsertSize");
                mReadWriter.write(",Duplicate,Secondary,Supplementary,NegStrand,ProperPair,TransExons,BestMatch,TransExonData");
                mReadWriter.newLine();
            }

            for(final ReadRecord read : reads)
            {
                mReadWriter.write(String.format("%s,%s,%s,%d,%d,%s,%d",
                        reads.size(), read.Id, read.Chromosome, read.PosStart, read.PosEnd, read.Cigar.toString(), read.fragmentInsertSize()));

                mReadWriter.write(String.format(",%s,%s,%s,%s,%s",
                    read.isDuplicate(), read.isSecondaryAlignment(), read.isSupplementaryAlignment(), read.isNegStrand(), read.isProperPair()));

                // log the transcript exons affected, and the highest matching transcript
                String transExonData = "";

                int transExonRefCount = 0;
                RegionMatchType highestTransMatchType = RegionMatchType.NONE;
                for(final Map.Entry<RegionMatchType,List<TransExonRef>> entry : read.getTransExonRefs().entrySet())
                {
                    for(final TransExonRef transExonRef : entry.getValue())
                    {
                        ++transExonRefCount;
                        transExonData = appendStr(transExonData, String.format("%s:%d:%s",
                                transExonRef.TransName, transExonRef.ExonRank, entry.getKey()), ';');

                        highestTransMatchType = entry.getKey();
                    }
                }

                mReadWriter.write(String.format(",%d,%s,%s",
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

}
