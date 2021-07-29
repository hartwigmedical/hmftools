package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;
import static com.hartwig.hmftools.neo.bind.BindData.DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class BlosumMapping
{
    private final int[][] mMappings;
    private boolean mIsValid;

    public static final String BLOSUM_FILE = "blosum_file";

    public BlosumMapping()
    {
        int aminoAcidCount = BindConstants.AMINO_ACIDS.size();
        mMappings = new int[aminoAcidCount][aminoAcidCount];
        mIsValid = false;
    }

    public boolean isValid() { return mIsValid; }

    public int map(final char aa1, final char aa2)
    {
        int aa1Index = aminoAcidIndex(aa1);
        int aa2Index = aminoAcidIndex(aa2);

        if(aa1Index == INVALID_AMINO_ACID || aa2Index == INVALID_AMINO_ACID)
            return 0;

        return mMappings[aa1Index][aa2Index];
    }

    public void load(final String filename)
    {
        if(filename == null || !Files.exists(Paths.get(filename)))
            return;

        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            String[] columns = lines.get(0).split(DELIM);
            lines.remove(0);

            if(columns.length != 21)
                return;

            Map<Integer,Integer> columnAaMap = Maps.newHashMap();

            for(int i = 1; i < columns.length; ++i)
            {
                char aa = columns[i].charAt(0);
                int aaIndex = aminoAcidIndex(aa);

                if(aaIndex == INVALID_AMINO_ACID)
                    return;

                columnAaMap.put(i, aaIndex);
            }

            for(String line : lines)
            {
                String[] items = line.split(DELIM);

                if(items.length != 21)
                    return;

                char aa1 = items[0].charAt(0);
                int aa1index = aminoAcidIndex(aa1);

                if(aa1index == INVALID_AMINO_ACID)
                    return;

                for(int i = 1; i < items.length; ++i)
                {
                    int correlation = Integer.parseInt(items[i]);
                    int aa2index = columnAaMap.get(i);

                    if(aa2index == INVALID_AMINO_ACID)
                        return;

                    mMappings[aa1index][aa2index] = correlation;
                }
            }

            NE_LOGGER.info("loaded blosum mappings from {}", filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load BLOSUM mapping file({}): {}", filename, e.toString());
        }

        mIsValid = true;
    }
}
