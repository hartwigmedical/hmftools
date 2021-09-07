package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DATA_TYPE_BIND_COUNTS;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;

import static org.apache.commons.math3.util.FastMath.log;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class FlankCounts
{
    // add in special placeholders for the start and end codons
    public static final char START_AMINO_ACID_ID = 'm';

    public static final double START_AMINO_ACID_FREQ = 0.00186249;
    public static final double STOP_AMINO_ACID_FREQ = START_AMINO_ACID_FREQ;

    public static final List<Character> FLANK_AMINO_ACIDS = Lists.newArrayList();

    public static final int FLANK_BASE_COUNT = 3;
    public static final int TOTAL_FLANK_BASE_COUNT = FLANK_BASE_COUNT * 2;

    public static final Map<Character,Integer> FLANK_AMINO_ACID_INDICES = Maps.newHashMap();

    static
    {
        FLANK_AMINO_ACIDS.addAll(AMINO_ACIDS);
        FLANK_AMINO_ACIDS.add(START_AMINO_ACID_ID);
        FLANK_AMINO_ACIDS.add(STOP_AMINO_ACID);

        for(int i = 0; i < FLANK_AMINO_ACIDS.size(); ++i)
        {
            FLANK_AMINO_ACID_INDICES.put(FLANK_AMINO_ACIDS.get(i), i);
        }
    }

    public static final int FLANK_AMINO_ACID_COUNT = FLANK_AMINO_ACIDS.size();

    public static int flankAminoAcidIndex(final char aminoAcid)
    {
        Integer index = FLANK_AMINO_ACID_INDICES.get(aminoAcid);
        return index != null ? index : INVALID_AMINO_ACID;
    }

    public static final int UP_1 = 2;
    public static final int DOWN_1 = 3;

    // by amino acid and position
    private final int[][] mBindCounts;

    private int mTotalBinds;
    private int mTotalWithFlanks;

    public FlankCounts()
    {
        mBindCounts = new int[FLANK_AMINO_ACID_COUNT][TOTAL_FLANK_BASE_COUNT];
        mTotalBinds = 0;
    }

    public final int[][] getBindCounts() { return mBindCounts; }

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
            ++mBindCounts[flankAminoAcidIndex(START_AMINO_ACID_ID)][UP_1];
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
            ++mBindCounts[flankAminoAcidIndex(STOP_AMINO_ACID)][DOWN_1];
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

        int aaIndex = flankAminoAcidIndex(aminoAcid);
        if(aaIndex == INVALID_AMINO_ACID)
            return;

        ++mBindCounts[aaIndex][flankPos];
    }

    public void writeData(final String filename, final FlankScores flankScores)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("DataType,AminoAcid");

            for(int i = FLANK_BASE_COUNT; i >= 1; --i)
            {
                writer.write(String.format(",U%d", i));
            }

            for(int i = 1; i <= FLANK_BASE_COUNT; ++i)
            {
                writer.write(String.format(",D%d", i));
            }

            writer.newLine();

            // write pos weights
            flankScores.writeData(writer);

            // write counts
            for(int aa = 0; aa < FLANK_AMINO_ACID_COUNT; ++aa)
            {
                char aminoAcid = FLANK_AMINO_ACIDS.get(aa);

                writer.write(String.format("%s,%c", DATA_TYPE_BIND_COUNTS, aminoAcid));

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
