package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class FlankCounts
{
    // add in special placeholders for the start and end codons
    private static final char START_AMINO_ACID = 'm';

    private static final List<Character> FLANK_AMINO_ACIDS = Lists.newArrayList();

    private static final int FLANK_BASE_COUNT = 3;
    private static final int TOTAL_FLANK_BASE_COUNT = FLANK_BASE_COUNT * 2;

    private static final Map<Character,Integer> FLANK_AMINO_ACID_INDICES = Maps.newHashMap();

    static
    {
        FLANK_AMINO_ACIDS.addAll(AMINO_ACIDS);
        FLANK_AMINO_ACIDS.add(START_AMINO_ACID);
        FLANK_AMINO_ACIDS.add(STOP_AMINO_ACID);

        for(int i = 0; i < FLANK_AMINO_ACIDS.size(); ++i)
        {
            FLANK_AMINO_ACID_INDICES.put(FLANK_AMINO_ACIDS.get(i), i);
        }
    }

    private static final int FLANK_AMINO_ACID_COUNT = FLANK_AMINO_ACIDS.size();

    private static int aminoAcidIndex(final char aminoAcid)
    {
        Integer index = FLANK_AMINO_ACID_INDICES.get(aminoAcid);
        return index != null ? index : INVALID_AMINO_ACID;
    }

    private static final int UP_1 = 2;
    private static final int DOWN_1 = 3;

    private final int[][] mBindCounts;
    private int mTotalBinds;
    private int mTotalWithFlanks;

    public FlankCounts()
    {
        mBindCounts = new int[FLANK_AMINO_ACID_COUNT][TOTAL_FLANK_BASE_COUNT];
        mTotalBinds = 0;
    }

    public void processBindData(final BindData bindData)
    {
        ++mTotalBinds;

        if(!bindData.hasFlanks())
            return;

        ++mTotalWithFlanks;

        // for up flank ACD and down flank EFG, the convention is to map them to 0-2 and 3-5
        if(!bindData.UpFlank.isEmpty())
        {
            int flankLength = bindData.UpFlank.length();

            if(flankLength >= 3)
            {
                addFlankBase(bindData.UpFlank.charAt(0), 0);
                addFlankBase(bindData.UpFlank.charAt(1), 1);
                addFlankBase(bindData.UpFlank.charAt(2), 2);
            }
            else if(flankLength >= 2)
            {
                addFlankBase(bindData.UpFlank.charAt(0), 1);
                addFlankBase(bindData.UpFlank.charAt(1), 2);
            }
            else
            {
                addFlankBase(bindData.UpFlank.charAt(0), 2);
            }
        }
        else
        {
            ++mBindCounts[aminoAcidIndex(START_AMINO_ACID)][UP_1];
        }

        if(!bindData.DownFlank.isEmpty())
        {
            int flankLength = bindData.DownFlank.length();

            if(flankLength >= 3)
            {
                addFlankBase(bindData.DownFlank.charAt(0), 3);
                addFlankBase(bindData.DownFlank.charAt(1), 4);
                addFlankBase(bindData.DownFlank.charAt(2), 5);
            }
            else if(flankLength >= 2)
            {
                addFlankBase(bindData.DownFlank.charAt(0), 3);
                addFlankBase(bindData.DownFlank.charAt(1), 4);
            }
            else
            {
                addFlankBase(bindData.DownFlank.charAt(0), 3);
            }
        }
        else
        {
            ++mBindCounts[aminoAcidIndex(STOP_AMINO_ACID)][DOWN_1];
        }
    }

    public void logStats()
    {
        NE_LOGGER.info("flanking {} from {} total binds", mTotalWithFlanks, mTotalBinds);
    }

    private void addFlankBase(char aminoAcid, int flankPos)
    {
        // explicitly ignore 21st
        if(aminoAcid == STOP_AMINO_ACID)
            return;

        int aaIndex = aminoAcidIndex(aminoAcid);
        if(aaIndex == INVALID_AMINO_ACID)
            return;

        ++mBindCounts[aaIndex][flankPos];
    }

    public void writeCounts(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("AminoAcid");

            for(int i = FLANK_BASE_COUNT; i >= 1; --i)
            {
                writer.write(String.format(",U%d", i));
            }

            for(int i = 1; i <= FLANK_BASE_COUNT; ++i)
            {
                writer.write(String.format(",D%d", i));
            }

            writer.newLine();

            for(int aa = 0; aa < FLANK_AMINO_ACID_COUNT; ++aa)
            {
                char aminoAcid = FLANK_AMINO_ACIDS.get(aa);

                writer.write(String.format("%c", aminoAcid));

                for(int pos = 0; pos < TOTAL_FLANK_BASE_COUNT; ++pos)
                {
                    writer.write(String.format(",%d", mBindCounts[aa][pos]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write flank counts data: {}", e.toString());
        }
    }
}
