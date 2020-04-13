package com.hartwig.hmftools.isofox.novel;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class FusionFinder
{
    private final IsofoxConfig mConfig;

    private final List<FusionReadData> mFusionReadData;

    public FusionFinder(final IsofoxConfig config)
    {
        mConfig = config;

        mFusionReadData = Lists.newArrayList();
    }

    public void findFusions(final Map<String,List<ChimericRead>> chimericReadsMap)
    {
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

}
