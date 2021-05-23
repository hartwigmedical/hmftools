package com.hartwig.hmftools.cobalt.diploid;

import static com.hartwig.hmftools.cobalt.CobaltConstants.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

class DiploidCount implements Comparable<DiploidCount>
{
    public final GenomePosition Position;
    
    private int mDiploid;
    private int mCount;

    private DiploidCount(@NotNull String line)
    {
        final String[] values = line.split(DELIMITER);
        Position = GenomePositions.create(values[0], Long.parseLong(values[1]));
        mDiploid = Integer.parseInt(values[2]);
        mCount = Integer.parseInt(values[3]);
    }

    public DiploidCount(final GenomePosition position, final int diploid, final int count)
    {
        Position = position;
        mDiploid = diploid;
        mCount = count;
    }

    void incrementTotal()
    {
        mCount++;
    }

    void incrementDiploid()
    {
        mDiploid++;
    }

    public int getDiploid()
    {
        return mDiploid;
    }

    public String chromosome()
    {
        return Position.chromosome();
    }

    public long position()
    {
        return Position.position();
    }

    public double proportionIsDiploid(int count)
    {
        return 1d * getDiploid() / count;
    }

    @NotNull
    public String toString()
    {
        return new StringJoiner(DELIMITER).add(Position.chromosome())
                .add(String.valueOf(Position.position()))
                .add(String.valueOf(mDiploid))
                .add(String.valueOf(mCount))
                .toString();
    }

    @Override
    public int compareTo(@NotNull final DiploidCount o)
    {
        return Position.compareTo(o.Position);
    }

    @NotNull
    public static Map<GenomePosition, DiploidCount> readDiploidCountAsMap(final String inputFile) throws IOException
    {
        return Files.readAllLines(new File(inputFile).toPath())
                .stream()
                .map(DiploidCount::new)
                .collect(Collectors.toMap(x -> x.Position, x -> x));
    }

    @NotNull
    public static List<DiploidCount> readDiploidCountAsList(final String inputFile) throws IOException
    {
        return Files.readAllLines(new File(inputFile).toPath()).stream().map(DiploidCount::new).collect(Collectors.toList());
    }
}
