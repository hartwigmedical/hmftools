package com.hartwig.hmftools.isofox.unmapped;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.util.StringJoiner;

public class BlatResult
{
    public final int Match;
    public final int Mismatch;
    public final int QGapCount;
    public final int QGapBases;
    public final int TGapCount;
    public final int TGapBases;
    public final byte Strand;
    public final String QName;
    public final int QSize;
    public final int QStart;
    public final int QEnd;
    public final int TSize;
    public final String TName;
    public final int TStart;
    public final int TEnd;
    public final int BlockCount;
    public final String BlockSizes;
    public final String QStarts;
    public final String TStarts;

    public BlatResult(final String data)
    {
        /*
        psLayout version 3
        0       1       2       3       4       5       6       7       8       9               10      11      12      13              14      15      16      17      18              19       20
        match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts
                match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
        ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        49      2       0       0       0       0       0       0       +       CPCT02010786T_65720606_12345    63      12      63      L1HS    6018    4081    4132    1       51,     12,     4081,
         */

        String[] values = data.split("\t", -1);

        Match = Integer.parseInt(values[0]);
        Mismatch = Integer.parseInt(values[1]);
        QGapCount = Integer.parseInt(values[4]);
        QGapBases = Integer.parseInt(values[5]);
        TGapCount = Integer.parseInt(values[6]);
        TGapBases = Integer.parseInt(values[7]);
        Strand = values[8].equals("+") ? POS_STRAND : NEG_STRAND;
        QName = values[9];
        QSize = Integer.parseInt(values[10]);
        QStart = Integer.parseInt(values[11]);
        QEnd = Integer.parseInt(values[12]);
        TName = values[13];
        TSize = Integer.parseInt(values[14]);
        TStart = Integer.parseInt(values[15]);
        TEnd = Integer.parseInt(values[16]);
        BlockCount = Integer.parseInt(values[17]);
        BlockSizes = values[18];
        QStarts = values[19];
        TStarts = values[20];
    }

    public static String csvHeader()
    {
        return "Match,Mismatch,Score,QGapCount,QGapBases,TGapCount,TGapBases,Strand,QStart,QEnd,TStart,TEnd,BlockCount,BlockSizes,QStarts,TStarts";

        /*
        StringJoiner sj = new StringJoiner(DELIMITER);
        sj.add()
        return sj.toString();
        */
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIMITER);
        sj.add(String.valueOf(Match));
        sj.add(String.valueOf(Mismatch));
        sj.add(String.valueOf(score()));
        sj.add(String.valueOf(QGapCount));
        sj.add(String.valueOf(QGapBases));
        sj.add(String.valueOf(TGapCount));
        sj.add(String.valueOf(TGapBases));
        sj.add(String.valueOf(Strand));
        sj.add(String.valueOf(QStart));
        sj.add(String.valueOf(QEnd));
        sj.add(String.valueOf(TStart));
        sj.add(String.valueOf(TEnd));
        sj.add(String.valueOf(BlockCount));
        sj.add(BlockSizes.replaceAll(",", ITEM_DELIM));
        sj.add(QStarts.replaceAll(",", ITEM_DELIM));
        sj.add(TStarts.replaceAll(",", ITEM_DELIM));
        return sj.toString();
    }

    public int score()
    {
        double qSizeMult = round(max(9 - floor(QSize/100), 1));
        int score = (int)round((1000 - (qSizeMult * Mismatch + QGapCount + TGapCount)) * min(Match/(double)QSize, 1));
        return score;
    }
}

