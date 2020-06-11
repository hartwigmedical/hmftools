package com.hartwig.hmftools.linx.annotators;

import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.jetbrains.annotations.Nullable;

public class ViralInsertAnnotator
{
    final private Map<String, String> mVirualHosts;

    public ViralInsertAnnotator()
    {
        mVirualHosts = Maps.newHashMap();
    }

    public void loadViralHostData(final String filename)
    {
        if(filename.isEmpty())
            return;

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine(); // skip header

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if(items.length < 2)
                {
                    LNX_LOGGER.debug("invalid data, entry {}", mVirualHosts.size());
                    break;
                }

                final String vhId = items[0];
                final String vhName = items[1];

                mVirualHosts.put(vhId, vhName);
            }

            LNX_LOGGER.debug("loaded {} viral host data", mVirualHosts.size());
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("failed to read viral host CSV file({})", filename);
        }
    }

    public static final int VH_ID = 0;
    public static final int VH_NAME = 1;

    @Nullable
    public String[] matchesViralInsert(final SvVarData svData)
    {
        if(mVirualHosts.isEmpty())
            return null;

        final String insertSeqAlignments = svData.getSvData().insertSequenceAlignments();
        if(insertSeqAlignments.isEmpty())
            return null;

        String[] vhData = insertSeqAlignments.split(":");
        if(vhData.length == 0)
            return null;

        final String vhId = vhData[0];
        final String vhName = mVirualHosts.get(vhId);

        if(vhName == null)
            return null;

        return new String[] {vhId, vhName};
    }


}
