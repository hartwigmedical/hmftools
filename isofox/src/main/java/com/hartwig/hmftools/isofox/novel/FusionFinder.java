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
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionFinder
{
    private final IsofoxConfig mConfig;

    private final List<FusionReadData> mFusionReadData;
    private BufferedWriter mReadWriter;

    public FusionFinder(final IsofoxConfig config)
    {
        mConfig = config;
        mFusionReadData = Lists.newArrayList();

        mReadWriter = null;
    }

    public void findFusions(final Map<String,List<ChimericRead>> chimericReadsMap)
    {
        // make a list of chimeric reads by their ID
        Map<String,List<ChimericRead>> readsByIdMap = Maps.newHashMap();

        for(final List<ChimericRead> chimericReads : chimericReadsMap.values())
        {
            for(final ChimericRead read : chimericReads)
            {
                List<ChimericRead> readsById = readsByIdMap.get(read.Id);

                if(readsById == null)
                {
                    readsById = Lists.newArrayList();
                    readsByIdMap.put(read.Id, readsById);
                }

                readsById.add(read);
            }
        }

        int unpairedReads = 0;

        for(Map.Entry<String,List<ChimericRead>> entry : readsByIdMap.entrySet())
        {
            final String readId = entry.getKey();

            final List<ChimericRead> reads = entry.getValue();

            if(reads.size() == 1)
            {
                ++unpairedReads;
                continue;
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

        ISF_LOGGER.info("chimeric fragments({}) unpaired({})", readsByIdMap.size(), unpairedReads);


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

    private ChimericRead findPairedRead(final ChimericRead read, final Map<String,List<ChimericRead>> chimericReadsMap)
    {
        final List<ChimericRead> chimericReads = chimericReadsMap.get(read.MateChromosome);

        if(chimericReads == null)
            return null;

        for(int index = 0; index < chimericReads.size(); ++index)
        {
            final ChimericRead otherRead = chimericReads.get(index);

            if(otherRead.Id.equals(read.Id))
            {
                chimericReads.remove(index);
                return otherRead;
            }
        }

        return null;
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

    private void writeReadData(final List<ChimericRead> reads)
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

            for(final ChimericRead read : reads)
            {
                mReadWriter.write(String.format("%s,%s,%s,%d,%d,%s,%d",
                        reads.size(), read.Id, read.Chromosome, read.PosStart, read.PosEnd, read.Cigar.toString(), read.InsertSize));

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
