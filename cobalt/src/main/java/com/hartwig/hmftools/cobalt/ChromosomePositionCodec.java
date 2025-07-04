package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import tech.tablesaw.api.LongColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.Table;

// tablesaw is not efficient at joining tables with chromosome, position columns.
// we therefore encode the chromosome, position into a long to help it.
// this is unfortunate, as pandas has no issue with the merge
public class ChromosomePositionCodec
{
    // chromosome number * 10_000_000_000 + position
    // we need to use long, longest chromosome is chr1, 242m bases
    private static final long CHROMOSOME_MULT = 10_000_000_000L;

    // to make sure conversion is consistent (1 vs chr1), we store the string we used for encoding
    // this is checked and retrieved for encoding. This will help catch problems with inconsistent
    // genome file versions
    private final Map<Long, String> mChromosomeNumStringMap;

    public ChromosomePositionCodec()
    {
        mChromosomeNumStringMap = new HashMap<>();
    }

    // this is used unfortunately due to tablesaw very slow merge when we merge by [chromosome, position]
    // pandas is able to complete in seconds
    public void addEncodedChrPosColumn(final Table table, boolean removeChromosomePosColumns)
    {
        LongColumn encodedChrPosCol = LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS);

        for(Row row : table)
        {
            String chr = row.getString(CobaltColumns.CHROMOSOME);

            int pos = row.getInt(CobaltColumns.POSITION);
            long chrPosIndex = encodeChromosomePosition(chr, pos);
            if(chrPosIndex >= 0)
            {
                encodedChrPosCol.append(chrPosIndex);
            }
            else
            {
                // add a missing value
                encodedChrPosCol.appendMissing();
            }
        }

        table.addColumns(encodedChrPosCol);

        if(removeChromosomePosColumns)
        {
            table.removeColumns(CobaltColumns.CHROMOSOME, CobaltColumns.POSITION);
        }
    }

    // add a chromosome position index that is the
    // chromosome number * 1_000_000_000 + position
    // we need to use long, longest chromosome is chr1, 242m bases
    public long encodeChromosomePosition(String chromosome, int position)
    {
        long chrNum;
        try
        {
            chrNum = HumanChromosome.fromString(chromosome).intValue();
        }
        catch (IllegalArgumentException e)
        {
            CB_LOGGER.fatal("unknown chromosome: {}", chromosome);
            throw new RuntimeException("unknown chromosome: " + chromosome);
        }

        // put it into the map, if fail raise exception, as it indicates chromosome name mismatch
        String chrNameCheck = mChromosomeNumStringMap.putIfAbsent(chrNum, chromosome);

        if(chrNameCheck != null && !chrNameCheck.equals(chromosome))
        {
            CB_LOGGER.fatal("inconsistent chromosome name: {} and {}", chrNameCheck, chromosome);
            throw new RuntimeException("inconsistent chromosome name");
        }

        long encoded = chrNum * CHROMOSOME_MULT + position;

        if(encoded < 0)
        {
            CB_LOGGER.fatal("chr: {} pos: {}, encoded: {} is negative",
                    chromosome, position, encoded);
            throw new RuntimeException("encoded chromosome position is negative");
        }

        return encoded;
    }

    // add a chromosome position index that is the
    // chromosome number * 10_000_000_000 + position
    // we need to use long, longest chromosome is chr1, 242m bases
    // i.e. chr3:22543005 would be encoded as 30_022_543_005
    public String decodeChromosome(long chromsomePositionCode)
    {
        long chrNum = chromsomePositionCode / CHROMOSOME_MULT;
        String chromosome = mChromosomeNumStringMap.get(chrNum);

        if(chromosome == null)
        {
            CB_LOGGER.error("unknown chromosome num: {}, encoded: {}",
                    chrNum, chromsomePositionCode);
        }

        return chromosome;
    }

    public int decodePosition(long chromsomePositionCode)
    {
        return (int)(chromsomePositionCode % CHROMOSOME_MULT);
    }
}
