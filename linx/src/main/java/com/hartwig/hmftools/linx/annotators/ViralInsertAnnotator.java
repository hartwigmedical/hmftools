package com.hartwig.hmftools.linx.annotators;

import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.jetbrains.annotations.Nullable;

public class ViralInsertAnnotator
{
    final private Map<String, String> mViralHosts;

    public ViralInsertAnnotator()
    {
        mViralHosts = Maps.newHashMap();
    }

    public void loadViralHostData(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());
            fileContents.remove(0); // assume a header

            for(final String line : fileContents)
            {
                String[] items = line.split(",");

                if(items.length < 2)
                {
                    LNX_LOGGER.debug("invalid data, entry {}", mViralHosts.size());
                    break;
                }

                final String vhId = items[0];
                final String vhName = items[1];

                mViralHosts.put(vhId, vhName);
            }

            LNX_LOGGER.debug("loaded {} viral host data", mViralHosts.size());
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
        if(mViralHosts.isEmpty())
            return null;

        final String insertSeqAlignments = svData.getSvData().insertSequenceAlignments();
        if(insertSeqAlignments.isEmpty())
            return null;

        String[] vhData = insertSeqAlignments.split(":");
        if(vhData.length == 0)
            return null;

        final String vhId = vhData[0];
        final String vhName = mViralHosts.get(vhId);

        if(vhName == null)
            return null;

        return new String[] {vhId, vhName};
    }


}
